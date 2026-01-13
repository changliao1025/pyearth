import numpy as np
from typing import Union, Tuple


def kde_support(
    data: np.ndarray,
    bandwidth: Union[float, str],
    grid_size: int,
    cut: float,
    bounds: np.ndarray,
) -> np.ndarray:
    """
    Establish support (range) for kernel density estimation.

    Calculates the evaluation points for kernel density estimation by extending
    the data range and constraining it within specified bounds. This creates
    a grid of points where the KDE will be evaluated.

    Parameters
    ----------
    data : np.ndarray
        Input data array (1D) for which to establish KDE support.
    bandwidth : Union[float, str]
        Bandwidth value used for extending the data range. Can be a numeric
        value or a string specifying bandwidth selection method.
    grid_size : int
        Number of equally spaced points in the support grid.
    cut : float
        Factor by which to extend the data range beyond min/max values.
        Extension amount = bandwidth * cut.
    bounds : np.ndarray
        Two-element array [min_bound, max_bound] constraining the support range.
        The calculated support will be clipped to these bounds.

    Returns
    -------
    np.ndarray
        1D array of equally spaced points defining the KDE support domain.
        Length equals grid_size.

    Notes
    -----
    The support range is calculated as:
    - Lower bound: max(data.min() - bandwidth * cut, bounds[0])
    - Upper bound: min(data.max() + bandwidth * cut, bounds[1])

    This ensures the KDE evaluation domain extends slightly beyond the data
    while respecting user-specified bounds to prevent excessive extrapolation.

    Examples
    --------
    Create support for univariate data:

    >>> import numpy as np
    >>> data = np.array([1, 2, 3, 4, 5])
    >>> support = kde_support(data, bandwidth=0.5, grid_size=100,
    ...                       cut=3.0, bounds=np.array([0, 10]))
    >>> len(support)
    100
    >>> support.min() >= 0 and support.max() <= 10
    True
    """
    # Calculate extended range beyond data limits
    data_min_extended = data.min() - bandwidth * cut
    data_max_extended = data.max() + bandwidth * cut

    # Constrain to user-specified bounds
    support_min = max(data_min_extended, bounds[0])
    support_max = min(data_max_extended, bounds[1])

    # Create equally spaced grid points
    return np.linspace(support_min, support_max, grid_size)


def scipy_bivariate_kde(
    x_data: np.ndarray,
    y_data: np.ndarray,
    bandwidth: Union[float, str],
    grid_size: int,
    cut: float,
    bounds: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute bivariate kernel density estimate using SciPy.

    Estimates the joint probability density function of two variables using
    Gaussian kernels. Creates a 2D grid of density values suitable for
    contour plotting or 3D surface visualization.

    Parameters
    ----------
    x_data : np.ndarray
        1D array of x-coordinates for the data points.
    y_data : np.ndarray
        1D array of y-coordinates for the data points. Must have same length as x_data.
    bandwidth : Union[float, str]
        Bandwidth for the kernel density estimation. Can be:
        - 'scott': Scott's rule of thumb
        - 'silverman': Silverman's rule of thumb
        - float: Fixed bandwidth value
    grid_size : int
        Number of points in each dimension of the evaluation grid.
    cut : float
        Factor by which to extend the data range for evaluation grid.
        Larger values provide smoother density estimates but may extrapolate beyond data.
    bounds : np.ndarray
        2x2 array [[x_min, x_max], [y_min, y_max]] constraining the evaluation domain.
        Format: np.array([[x_bounds], [y_bounds]]) where each subarray is [min, max].

    Returns
    -------
    Tuple[np.ndarray, np.ndarray, np.ndarray]
        - xx : np.ndarray
            2D meshgrid of x-coordinates (grid_size x grid_size)
        - yy : np.ndarray
            2D meshgrid of y-coordinates (grid_size x grid_size)
        - z : np.ndarray
            2D array of density values at each (x,y) grid point

    Raises
    ------
    ImportError
        If scipy is not installed.
    ValueError
        If bandwidth parameter cannot be handled by scipy backend, or if
        different bandwidths are specified for each dimension.

    Notes
    -----
    This function uses scipy.stats.gaussian_kde for kernel density estimation.
    The evaluation grid is created using meshgrid from the support ranges
    calculated for each dimension.

    Bandwidth handling:
    - String values ('scott', 'silverman') use scipy's built-in bandwidth selection
    - Scalar values apply the same bandwidth to both dimensions
    - Different bandwidths per dimension are not supported (use statsmodels backend)

    The density values in z are normalized such that their integral over the
    domain approximates 1 (probability density function).

    Examples
    --------
    Basic bivariate KDE:

    >>> import numpy as np
    >>> x = np.random.normal(0, 1, 100)
    >>> y = np.random.normal(0, 1, 100)
    >>> bounds = np.array([[x.min(), x.max()], [y.min(), y.max()]])
    >>> xx, yy, z = scipy_bivariate_kde(x, y, 'scott', 50, 3.0, bounds)
    >>> xx.shape, yy.shape, z.shape
    ((50, 50), (50, 50), (50, 50))

    Custom bandwidth:

    >>> xx, yy, z = scipy_bivariate_kde(x, y, bandwidth=0.1, grid_size=100,
    ...                                cut=2.0, bounds=bounds)

    See Also
    --------
    kde_support : Calculate support range for univariate KDE
    scipy.stats.gaussian_kde : Underlying SciPy KDE implementation
    """
    try:
        import scipy
    except ImportError as e:
        raise ImportError(
            "The package 'scipy' is required for this function to run."
        ) from e

    # Combine x and y data into coordinate pairs
    data = np.c_[x_data, y_data]

    # Create KDE object
    kde = scipy.stats.gaussian_kde(data.T, bw_method=bandwidth)

    # Calculate data standard deviations for bandwidth scaling
    data_std = data.std(axis=0, ddof=1)

    # Handle different bandwidth specifications
    if isinstance(bandwidth, str):
        # Correct common typo and get bandwidth factors
        bw_method = "scotts" if bandwidth == "scott" else bandwidth
        bw_x = getattr(kde, f"{bw_method}_factor")() * data_std[0]
        bw_y = getattr(kde, f"{bw_method}_factor")() * data_std[1]
    elif np.isscalar(bandwidth):
        # Same bandwidth for both dimensions
        bw_x = bw_y = bandwidth
    else:
        # Different bandwidths not supported by scipy backend
        msg = (
            "Cannot specify a different bandwidth for each dimension "
            "with the scipy backend. You should install statsmodels."
        )
        raise ValueError(msg)

    # Create evaluation grids for x and y dimensions
    x_support = kde_support(data[:, 0], bw_x, grid_size, cut, bounds[0])
    y_support = kde_support(data[:, 1], bw_y, grid_size, cut, bounds[1])

    # Create 2D meshgrid for evaluation points
    xx, yy = np.meshgrid(x_support, y_support)

    # Evaluate KDE at all grid points
    z = kde([xx.ravel(), yy.ravel()]).reshape(xx.shape)

    return xx, yy, z
