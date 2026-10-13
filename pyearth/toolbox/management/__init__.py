"""Dataset management tools for raster and vector geospatial files."""

from pyearth.toolbox.management.raster.merge_raster_files import (
    merge_raster_files,
)
from pyearth.toolbox.management.raster.reproject import reproject_raster
from pyearth.toolbox.management.raster.resample import resample_raster
from pyearth.toolbox.management.vector.fields import (
    add_field_to_vector_file,
    get_field_and_value,
    get_field_value,
)
from pyearth.toolbox.management.vector.merge_features import (
    merge_features,
)
from pyearth.toolbox.management.vector.merge_vector_files import (
    merge_vector_files,
)
from pyearth.toolbox.management.vector.polygon_calculator import (
    polygon_difference,
)
from pyearth.toolbox.management.vector.remove_small_polygon import (
    remove_small_polygon,
)
from pyearth.toolbox.management.vector.reproject import (
    reproject_vector,
)

__all__ = [
    "merge_raster_files",
    "reproject_raster",
    "resample_raster",
    "merge_vector_files",
    "merge_features",
    "reproject_vector",
    "remove_small_polygon",
    "get_field_value",
    "get_field_and_value",
    "add_field_to_vector_file",
    "polygon_difference",
]
