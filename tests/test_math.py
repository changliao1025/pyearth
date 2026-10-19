import numpy as np
import pytest

from pyearth.toolbox.math import (
    gap_fill_by_window,
    kde_support,
    remap,
    scipy_bivariate_kde,
    search_neighbors,
)


class TestRemap:
    def test_basic_remap(self):
        assert remap(5, 0, 10, 0, 100) == 50.0
        assert remap(0, 0, 10, 10, 20) == 10.0
        assert remap(10, 0, 10, 10, 20) == 20.0

    def test_reversed_ranges(self):
        # Reversed output
        assert remap(5, 0, 10, 100, 0) == 50.0
        # Reversed input
        assert remap(5, 10, 0, 0, 100) == 50.0
        # Both reversed
        assert remap(5, 10, 0, 100, 0) == 50.0

    def test_edge_cases(self):
        # Zero input range returns None
        assert remap(5, 10, 10, 0, 100) is None
        # Zero output range returns constant nMin
        assert remap(5, 0, 10, 50, 50) == 50.0


class TestGapFillByWindow:
    def test_search_neighbors(self):
        arr = np.array([10.0, 20.0, 30.0, 40.0, 50.0])
        # Center index 2 with window 1 -> indices 1, 2, 3
        nb = search_neighbors(2, arr, 1)
        assert np.array_equal(nb, [20.0, 30.0, 40.0])

        # Boundary index 0 with window 2 -> indices 0, 1, 2
        nb_edge = search_neighbors(0, arr, 2)
        assert np.array_equal(nb_edge, [10.0, 20.0, 30.0])

    def test_search_neighbors_validation(self):
        arr = np.array([1.0, 2.0])
        with pytest.raises(ValueError):
            search_neighbors(5, arr, 1)
        with pytest.raises(ValueError):
            search_neighbors(0, arr, -1)
        with pytest.raises(ValueError):
            search_neighbors(0, np.ones((2, 2)), 1)

    def test_gap_fill_by_window(self):
        arr = np.array([1.0, np.nan, 3.0, np.nan, 5.0])
        filled = gap_fill_by_window(arr, iWindow_size_in=1, min_good_neighbors=1)
        assert np.isclose(filled[1], 2.0)
        assert np.isclose(filled[3], 4.0)

    def test_gap_fill_insufficient_neighbors(self):
        arr = np.array([np.nan, 2.0, 3.0])
        filled = gap_fill_by_window(arr, iWindow_size_in=1, min_good_neighbors=2)
        # First gap at index 0 only has neighbor index 1 (one good neighbor), so remains NaN
        assert np.isnan(filled[0])


class TestBivariateKDE:
    def test_kde_support(self):
        data = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        bounds = np.array([0.0, 10.0])
        support = kde_support(data, bandwidth=0.5, grid_size=50, cut=2.0, bounds=bounds)
        assert len(support) == 50
        assert support.min() >= 0.0
        assert support.max() <= 10.0

    def test_scipy_bivariate_kde(self):
        np.random.seed(42)
        x = np.random.normal(0, 1, 60)
        y = np.random.normal(0, 1, 60)
        bounds = np.array([[x.min() - 1, x.max() + 1], [y.min() - 1, y.max() + 1]])
        xx, yy, z = scipy_bivariate_kde(
            x, y, bandwidth="scott", grid_size=20, cut=1.0, bounds=bounds
        )
        assert xx.shape == (20, 20)
        assert yy.shape == (20, 20)
        assert z.shape == (20, 20)
        assert np.all(z >= 0.0)
