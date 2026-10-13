"""Data processing toolbox: time series conversion, outlier removal, percentiles, lists, vector geometry copy, ocean masks."""

from pyearth.toolbox.data.cgpercentiles import cgpercentiles
from pyearth.toolbox.data.convert_time_series_daily_to_monthly import (
    convert_time_series_daily_to_monthly,
)
from pyearth.toolbox.data.list.check_if_duplicates import check_if_duplicates
from pyearth.toolbox.data.list.list_alg import list_max, list_min
from pyearth.toolbox.data.ocean.define_land_ocean_mask import (
    create_land_ocean_vector_mask_naturalearth,
)
from pyearth.toolbox.data.remove_outliers import remove_outliers
from pyearth.toolbox.data.vector.copy_geometry_without_attributes import (
    copy_geometry_without_attributes,
)

__all__ = [
    "cgpercentiles",
    "convert_time_series_daily_to_monthly",
    "remove_outliers",
    "check_if_duplicates",
    "list_max",
    "list_min",
    "copy_geometry_without_attributes",
    "create_land_ocean_vector_mask_naturalearth",
]
