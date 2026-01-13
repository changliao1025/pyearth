"""
Create a rectangular latitude/longitude based mesh using GIS operations.

This module creates mesh cells representing edges instead of centers, using GDAL API
for most GIS operations. The mesh is defined by longitude left, latitude bottom,
number of rows/columns, and resolution.
"""

import os
from typing import List, Optional, Tuple, Union
import numpy as np
from osgeo import ogr, osr
from pyearth.gis.gdal.gdal_vector_format_support import get_vector_driver_from_filename
from pyearth.toolbox.mesh.algorithm.convert_coordinates import (
    convert_gcs_coordinates_to_meshcell,
)

# Constants
DEFAULT_EPSG_CODE = 4326  # WGS84 lat/lon
DEFAULT_PRECISION = 2
DEFAULT_FIELD_WIDTH = 20
COORDINATE_FILL_VALUE = -9999.0


def _validate_mesh_parameters(
    longitude_left: float,
    latitude_bottom: float,
    resolution_degrees: float,
    num_columns: int,
    num_rows: int,
    output_filename: str,
) -> None:
    """
    Validate input parameters for mesh creation.

    Args:
        longitude_left: Left boundary longitude
        latitude_bottom: Bottom boundary latitude
        resolution_degrees: Cell resolution in degrees
        num_columns: Number of columns
        num_rows: Number of rows
        output_filename: Output file path

    Raises:
        ValueError: If any parameter is invalid
        FileNotFoundError: If output directory doesn't exist
    """
    if not (-180.0 <= longitude_left <= 180.0):
        raise ValueError(
            f"Longitude must be between -180 and 180 degrees, got {longitude_left}"
        )

    if not (-90.0 <= latitude_bottom <= 90.0):
        raise ValueError(
            f"Latitude must be between -90 and 90 degrees, got {latitude_bottom}"
        )

    if resolution_degrees <= 0:
        raise ValueError(f"Resolution must be positive, got {resolution_degrees}")

    if num_columns <= 0 or num_rows <= 0:
        raise ValueError(
            f"Columns and rows must be positive integers, got columns={num_columns}, rows={num_rows}"
        )

    if not isinstance(output_filename, str) or not output_filename.strip():
        raise ValueError("Output filename must be a non-empty string")

    # Check if output directory exists
    output_dir = os.path.dirname(output_filename)
    if output_dir and not os.path.exists(output_dir):
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")


def _create_boundary_geometry(
    longitude_left: float,
    latitude_bottom: float,
    num_columns: int,
    num_rows: int,
    resolution_degrees: float,
) -> ogr.Geometry:
    """
    Create a boundary geometry for the mesh using bounding box.

    Args:
        longitude_left: Left boundary longitude
        latitude_bottom: Bottom boundary latitude
        num_columns: Number of columns
        num_rows: Number of rows
        resolution_degrees: Cell resolution in degrees

    Returns:
        ogr.Geometry: Polygon geometry representing the boundary
    """
    boundary = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)

    # Calculate bounding box coordinates
    x1, y1 = longitude_left, latitude_bottom
    x2, y2 = x1 + num_columns * resolution_degrees, y1
    x3, y3 = x2, y1 + num_rows * resolution_degrees
    x4, y4 = x1, y3

    # Add points in counter-clockwise order
    ring.AddPoint(x1, y1)
    ring.AddPoint(x2, y2)
    ring.AddPoint(x3, y3)
    ring.AddPoint(x4, y4)
    ring.AddPoint(x1, y1)  # Close the ring
    ring.CloseRings()

    boundary.AddGeometry(ring)
    return boundary


def create_latlon_mesh(
    dLongitude_left_in: float,
    dLatitude_bot_in: float,
    dResolution_degree_in: float,
    ncolumn_in: int,
    nrow_in: int,
    sFilename_output_in: str,
    pBoundary_in: Optional[str] = None,
) -> List:
    """
    Create a rectangular latitude/longitude based mesh with cell topology and neighbor relationships.

    This function generates a structured mesh grid using geographic coordinates (latitude/longitude)
    and creates mesh cells with their connectivity information. Each cell includes area calculations,
    neighbor relationships, and distance calculations between neighboring cells.

    Args:
        dLongitude_left_in (float): Left (western) boundary longitude in degrees (-180 to 180).
        dLatitude_bot_in (float): Bottom (southern) boundary latitude in degrees (-90 to 90).
        dResolution_degree_in (float): Cell resolution in degrees (must be positive).
        ncolumn_in (int): Number of columns in the mesh grid (must be positive integer).
        nrow_in (int): Number of rows in the mesh grid (must be positive integer).
        sFilename_output_in (str): Path to output GeoJSON file for mesh visualization.
        pBoundary_in (Optional[str]): Well-Known Text (WKT) string defining the boundary polygon.
                                    If None, a bounding box boundary will be created automatically.

    Returns:
        List: A list of mesh cell objects with the following properties:
            - lCellID: Unique cell identifier
            - pVertex_center: Cell center coordinates
            - aNeighbor: List of neighboring cell IDs
            - nNeighbor: Number of neighbors
            - aNeighbor_distance: List of distances to neighboring cells
            - Cell area and edge length calculations

    Raises:
        ValueError: If input parameters are invalid (negative values, invalid coordinates).
        FileNotFoundError: If output directory doesn't exist.
        RuntimeError: If GDAL operations fail or mesh generation encounters errors.

    Example:
        >>> mesh_cells = create_latlon_mesh(
        ...     dLongitude_left_in=-180.0,
        ...     dLatitude_bot_in=-90.0,
        ...     dResolution_degree_in=1.0,
        ...     ncolumn_in=360,
        ...     nrow_in=180,
        ...     sFilename_output_in="global_mesh.geojson"
        ... )
        >>> print(f"Created {len(mesh_cells)} mesh cells")

    Note:
        - Mesh cells are created using counter-clockwise vertex ordering
        - Cell indexing starts from lower-left corner (1-based)
        - Neighbor relationships are calculated using 8-connectivity
        - Output file is saved in GeoJSON format with WGS84 projection
        - For geometry objects, WKT strings are used to avoid crashes when datasets are closed
          (see https://gdal.org/api/python_gotchas.html)
    """
    # Input validation
    if not (-180.0 <= dLongitude_left_in <= 180.0):
        raise ValueError(
            f"Longitude must be between -180 and 180 degrees, got {dLongitude_left_in}"
        )

    if not (-90.0 <= dLatitude_bot_in <= 90.0):
        raise ValueError(
            f"Latitude must be between -90 and 90 degrees, got {dLatitude_bot_in}"
        )

    if dResolution_degree_in <= 0:
        raise ValueError(f"Resolution must be positive, got {dResolution_degree_in}")

    if ncolumn_in <= 0 or nrow_in <= 0:
        raise ValueError(
            f"Columns and rows must be positive integers, got columns={ncolumn_in}, rows={nrow_in}"
        )

    if not isinstance(sFilename_output_in, str) or not sFilename_output_in.strip():
        raise ValueError("Output filename must be a non-empty string")

    # Check if output directory exists
    output_dir = os.path.dirname(sFilename_output_in)
    if output_dir and not os.path.exists(output_dir):
        raise FileNotFoundError(f"Output directory does not exist: {output_dir}")
    # Handle boundary geometry
    # For geometry objects, WKT strings are used to avoid crashes when datasets are closed
    # Reference: https://gdal.org/api/python_gotchas.html
    if pBoundary_in is None:
        print("Creating mesh with automatic bounding box boundary")
        try:
            pBoundary = _create_boundary_geometry(
                dLongitude_left_in,
                dLatitude_bot_in,
                ncolumn_in,
                nrow_in,
                dResolution_degree_in,
            )
        except Exception as e:
            raise RuntimeError(f"Failed to create boundary geometry: {e}")
    else:
        try:
            pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)
            if pBoundary is None:
                raise ValueError(f"Invalid WKT string: {pBoundary_in}")
        except Exception as e:
            raise ValueError(f"Failed to create geometry from WKT: {e}")

    # Setup output file
    if os.path.exists(sFilename_output_in):
        try:
            os.remove(sFilename_output_in)
        except OSError as e:
            raise RuntimeError(f"Failed to remove existing output file: {e}")

    # Create output dataset and layer
    try:
        pDriver = get_vector_driver_from_filename(sFilename_output_in)
        if pDriver is None:
            raise RuntimeError("GeoJSON driver not available")

        pDataset = pDriver.CreateDataSource(sFilename_output_in)
        if pDataset is None:
            raise RuntimeError(
                f"Failed to create output dataset: {sFilename_output_in}"
            )

        # Setup spatial reference system (WGS84)
        pSpatial_reference_gcs = osr.SpatialReference()
        pSpatial_reference_gcs.ImportFromEPSG(DEFAULT_EPSG_CODE)
        pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

        # Create layer with attributes
        pLayer = pDataset.CreateLayer("cell", pSpatial_reference_gcs, ogr.wkbPolygon)
        if pLayer is None:
            raise RuntimeError("Failed to create layer")

        # Define attribute fields
        pLayer.CreateField(
            ogr.FieldDefn("cellid", ogr.OFTInteger64)
        )  # 64-bit for high resolution
        pLayer.CreateField(ogr.FieldDefn("longitude", ogr.OFTReal))
        pLayer.CreateField(ogr.FieldDefn("latitude", ogr.OFTReal))

        pArea_field = ogr.FieldDefn("area", ogr.OFTReal)
        pArea_field.SetWidth(DEFAULT_FIELD_WIDTH)
        pArea_field.SetPrecision(DEFAULT_PRECISION)
        pLayer.CreateField(pArea_field)

        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)

    except Exception as e:
        raise RuntimeError(f"Failed to setup output layer: {e}")
    xleft = dLongitude_left_in
    xspacing = dResolution_degree_in
    ybottom = dLatitude_bot_in
    yspacing = dResolution_degree_in
    aLatlon = list()
    aLatlon_dict = dict()
    lCellIndex = 0

    def add_cell_into_list(
        aList, lCellID, iRow, iColumn, dLongitude_center, dLatitude_center, aCoords
    ):
        pLatlon = convert_gcs_coordinates_to_meshcell(
            2, dLongitude_center, dLatitude_center, aCoords
        )
        pLatlon.lCellID = lCellID
        dArea = pLatlon.calculate_polygon_area()
        pLatlon.calculate_line_length()
        # build topoloy
        aNeighbor = list()
        aNeighbor_distance = list()
        # lCellID_center = lCellID
        # counter-clock wise direction to add the neighbor
        if iRow > 1:  # under
            iRow_dummy = iRow - 1
            if iColumn > 1:
                iColumn_dummy = iColumn - 1
                lCellID2 = (
                    iRow_dummy - 1
                ) * ncolumn_in + iColumn_dummy  # lCellID0 - nrow_in
                aNeighbor.append(lCellID2)

            lCellID0 = (iRow_dummy - 1) * ncolumn_in + iColumn
            aNeighbor.append(lCellID0)
        if iColumn < ncolumn_in:  # right
            iColumn_dummy = iColumn + 1
            if iRow > 1:
                iRow_dummy = iRow - 1
                lCellID7 = (iRow_dummy - 1) * ncolumn_in + iColumn_dummy  # lCellID5 -1
                aNeighbor.append(lCellID7)
            lCellID5 = (
                iRow - 1
            ) * ncolumn_in + iColumn_dummy  # nrow_in * iColumn + iRow
            aNeighbor.append(lCellID5)
        if iRow < nrow_in:  # top
            iRow_dummy = iRow + 1
            if iColumn < ncolumn_in:
                iColumn_dummy = iColumn + 1
                lCellID6 = (
                    iRow_dummy - 1
                ) * ncolumn_in + iColumn_dummy  # lCellID3 + nrow_in
                aNeighbor.append(lCellID6)
            lCellID3 = (iRow_dummy - 1) * ncolumn_in + iColumn  # lCellID_center + 1
            aNeighbor.append(lCellID3)

        if iColumn > 1:  # left
            iColumn_dummy = iColumn - 1
            if iRow < nrow_in:
                iRow_dummy = iRow + 1
                lCellID4 = (iRow_dummy - 1) * ncolumn_in + iColumn_dummy  # lCellID1 + 1
                aNeighbor.append(lCellID4)
            lCellID1 = (
                iRow - 1
            ) * ncolumn_in + iColumn_dummy  # nrow_in * (iColumn-2) + iRow
            aNeighbor.append(lCellID1)

        pLatlon.aNeighbor = aNeighbor
        pLatlon.nNeighbor = len(aNeighbor)
        pLatlon.aNeighbor_land = aNeighbor
        pLatlon.nNeighbor_land = pLatlon.nNeighbor
        aList.append(pLatlon)

        return aList, dArea

    # change the order because mpas uses counter-clock wise to store the vertices
    # we will also start from the lower-left corner, and then go to the right and then go up
    # so the final index will be like this
    # 3 4
    # 1 2
    # lCellID = 1
    # .........
    # (x4,y4)-----(x3,y3)
    #   |           |
    # (x1,y1)-----(x2,y2)
    # ...............

    try:
        for iRow in range(1, nrow_in + 1):
            for iColumn in range(1, ncolumn_in + 1):
                lCellID = (iRow - 1) * ncolumn_in + iColumn
                # define a polygon here
                x1 = xleft + ((iColumn - 1) * xspacing)
                y1 = ybottom + ((iRow - 1) * yspacing)

                x2 = xleft + ((iColumn) * xspacing)
                y2 = ybottom + ((iRow - 1) * yspacing)

                x3 = xleft + ((iColumn) * xspacing)
                y3 = ybottom + ((iRow) * yspacing)

                x4 = xleft + ((iColumn - 1) * xspacing)
                y4 = ybottom + ((iRow) * yspacing)

                coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]

                ring = ogr.Geometry(ogr.wkbLinearRing)
                try:
                    for x, y in coordinates:
                        ring.AddPoint(x, y)
                except Exception as e:
                    print(f"Error adding points to ring geometry: {e}")
                    raise

                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)
                aCoords_gcs = np.full((5, 2), -9999.0, dtype=float)
                try:
                    for i, (x, y) in enumerate(coordinates):
                        aCoords_gcs[i, 0] = x
                        aCoords_gcs[i, 1] = y
                except Exception as e:
                    print(f"Error setting coordinate arrays: {e}")
                    raise

                dLongitude_center = np.mean(aCoords_gcs[0:4, 0])
                dLatitude_center = np.mean(aCoords_gcs[0:4, 1])

                iFlag = False
                if pPolygon.Within(pBoundary):
                    iFlag = True
                else:
                    # then check intersection
                    if pPolygon.Intersects(pBoundary):
                        iFlag = True
                    else:
                        pass

                if iFlag == True:
                    aLatlon, dArea = add_cell_into_list(
                        aLatlon,
                        lCellID,
                        iRow,
                        iColumn,
                        dLongitude_center,
                        dLatitude_center,
                        aCoords_gcs,
                    )
                    # save feature
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", lCellID)
                    pFeature.SetField("longitude", dLongitude_center)
                    pFeature.SetField("latitude", dLatitude_center)
                    pFeature.SetField("area", dArea)
                    pLayer.CreateFeature(pFeature)
                    # add to dictionary
                    aLatlon_dict[lCellID] = lCellIndex
                    lCellIndex = lCellIndex + 1
                    pass
    except Exception as e:
        print(f"Error processing mesh grid rows and columns: {e}")
        raise

    pDataset = pLayer = pFeature = None

    aLatlon_out = list()
    try:
        for pCell in aLatlon:
            aNeighbor = pCell.aNeighbor
            aNeighbor_land_update = list()
            try:
                for lNeighbor in aNeighbor:
                    if lNeighbor in aLatlon_dict:
                        aNeighbor_land_update.append(lNeighbor)
            except Exception as e:
                print(f"Error processing neighbor list: {e}")
                raise
            # for latlon, there is no ocean concept
            pCell.aNeighbor = aNeighbor_land_update
            pCell.nNeighbor = len(aNeighbor_land_update)
            pCell.aNeighbor_land = aNeighbor_land_update
            pCell.nNeighbor_land = len(aNeighbor_land_update)
            pCell.nNeighbor_ocean = pCell.nPoint - pCell.nNeighbor_land
            aLatlon_out.append(pCell)
    except Exception as e:
        print(f"Error updating cell neighbors: {e}")
        raise

    # calculate neighbor distance
    try:
        for pLatlon in aLatlon_out:
            aNeighbor = pLatlon.aNeighbor
            pLatlon.aNeighbor_distance = list()
            try:
                for lCellID1 in aNeighbor:
                    # use dictionary to get index
                    lIndex = aLatlon_dict[lCellID1]
                    pLatlon1 = aLatlon_out[lIndex]
                    dDistance = pLatlon.pPoint_center.calculate_distance(
                        pLatlon1.pPoint_center
                    )
                    pLatlon.aNeighbor_distance.append(dDistance)
            except Exception as e:
                print(f"Error calculating neighbor distances: {e}")
                raise
    except Exception as e:
        print(f"Error processing neighbor distance calculation: {e}")
        raise

    return aLatlon_out
