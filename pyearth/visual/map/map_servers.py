import os
import numpy as np

from io import BytesIO
import math
from osgeo import osr
import cartopy.io.img_tiles as cimgt
from pyearth.gis.spatialref.convert_between_degree_and_meter import degree_to_meter

def lonlat_to_tile(lon, lat, zoom):
    """
    Convert longitude and latitude to tile indices at a given zoom level.
    """
    lat_rad = math.radians(lat)
    n = 2.0 ** zoom
    x_tile = int((lon + 180.0) / 360.0 * n)
    y_tile = int((1.0 - math.log(math.tan(lat_rad) + (1 / math.cos(lat_rad))) / math.pi) / 2.0 * n)
    return x_tile, y_tile

def extent_to_tile_indices(minx, miny, maxx, maxy, zoom):
    """
    Convert extent (minx, miny, maxx, maxy) to tile indices at a given zoom level.
    """
    x_min, y_max = lonlat_to_tile(minx, miny, zoom)
    x_max, y_min = lonlat_to_tile(maxx, maxy, zoom)
    return x_min, y_min, x_max, y_max

def combine_tiles(tiles, tile_size):
    from PIL import Image
    rows = len(tiles)
    cols = len(tiles[0])
    combined_img = Image.new('RGBA', (cols * tile_size, rows * tile_size))  # Use RGBA mode for transparency

    for row in range(rows):
        for col in range(cols):
            combined_img.paste(tiles[row][col], (col * tile_size, row * tile_size))

    return combined_img

class StadiaStamen(cimgt.Stamen):
    def _image_url(self, tile):
        x,y,z = tile
        #get sStadia_API_KEY from the environment variable
        sStadia_API_KEY = os.environ.get('STADIA_API_KEY')
        url = f"https://tiles.stadiamaps.com/tiles/stamen_terrain/{z}/{x}/{y}.png?api_key={sStadia_API_KEY}"
        return url

# Create a custom tile source for the Esri World Terrain Base
class EsriRelief(cimgt.GoogleTiles):
    def _image_url(self, tile):
        x, y, z = tile
        return f'https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg'

class EsriTerrain(cimgt.GoogleTiles):
    def _image_url(self, tile):
        x, y, z = tile
        return f'https://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}'

# Create a custom tile source for the Esri Hydro Reference Overlay
class EsriHydro(cimgt.GoogleTiles):
    def _image_url(self, tile):
        x, y, z = tile
        return f'https://tiles.arcgis.com/tiles/P3ePLMYs2RVChkJx/arcgis/rest/services/Esri_Hydro_Reference_Overlay/MapServer/tile/{z}/{y}/{x}'

def fetch_stadia_tile(z, x, y):
    import requests
    from PIL import Image
        #get sStadia_API_KEY from the environment variable
    sStadia_API_KEY = os.environ.get('STADIA_API_KEY')
    url = f"https://tiles.stadiamaps.com/tiles/stamen_terrain/{z}/{x}/{y}.png?api_key={sStadia_API_KEY}"
    response = requests.get(url)
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content))
        print("stadia tile size:",img.size)
        #save the image as a png file using x, y, z
        #img.save(os.path.join(sWorkspace_png, f'{x}_{y}_{z}.png'))
        return img
    else:
        raise Exception(f"Failed to fetch tile: {response.status_code}")

def fetch_esri_terrain_tile(z, x, y):
    import requests
    from PIL import Image
    url = f"https://server.arcgisonline.com/ArcGIS/rest/services/World_Terrain_Base/MapServer/tile/{z}/{y}/{x}"
    response = requests.get(url)
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content))
        #print tile image size
        print("esri terrain tile size:",img.size)
        #save the image as a png file using x, y, z
        #img.save(os.path.join(sWorkspace_png, f'{x}_{y}_{z}.png'))
        return img
    else:
        raise Exception(f"Failed to fetch tile: {response.status_code}")

def fetch_esri_relif_tile(z, x, y):
    import requests
    from PIL import Image
    url = f"https://server.arcgisonline.com/ArcGIS/rest/services/World_Shaded_Relief/MapServer/tile/{z}/{y}/{x}.jpg"
    response = requests.get(url)
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content))
        print("esri relief tile size:",img.size)
        return img
    else:
        raise Exception(f"Failed to fetch tile: {response.status_code}")

def fetch_esri_hydro_tile(z, x, y):
    import requests
    from PIL import Image
    url = f"https://tiles.arcgis.com/tiles/P3ePLMYs2RVChkJx/arcgis/rest/services/Esri_Hydro_Reference_Overlay/MapServer/tile/{z}/{y}/{x}"
    response = requests.get(url)
    if response.status_code == 200:
        img = Image.open(BytesIO(response.content)).convert('RGBA')  # Convert to RGBA mode
        print("esri hydro tile size:",img.size)
        datas = img.getdata()
        new_data = []
        for item in datas:
            # Change all black (also shades of black)
            # to transparent
            if item[0] == 0 and item[1] == 0 and item[2] == 0:
                new_data.append((0, 0, 0, 0))
            else:
                new_data.append(item)
        img.putdata(new_data)
        return img
    else:
        raise Exception(f"Failed to fetch tile: {response.status_code}")

def Stadia_terrain_images(aExtent, zoom_level):
    minx, maxx, miny,  maxy = aExtent
    x_min, y_min, x_max, y_max = extent_to_tile_indices(minx, miny, maxx, maxy, zoom_level)
    tile_size = 256
    tiles = []
    for y in range(y_min, y_max + 1):
        row = []
        for x in range(x_min, x_max + 1):
            tile = fetch_stadia_tile(zoom_level, x, y)
            row.append(tile)
        tiles.append(row)

    # Combine the tiles into a single image
    combined_img = combine_tiles(tiles, tile_size)

    # Convert the image to a NumPy array
    img_array = np.array(combined_img)
    return img_array

def Esri_terrain_images(aExtent, zoom_level):
    minx, maxx, miny, maxy = aExtent
    x_min, y_min, x_max, y_max = extent_to_tile_indices(minx, miny, maxx, maxy, zoom_level)
    tile_size = 256
    tiles = []
    for y in range(y_min, y_max + 1):
        row = []
        for x in range(x_min, x_max + 1):
            tile = fetch_esri_terrain_tile(zoom_level, x, y)
            row.append(tile)
        tiles.append(row)

    # Combine the tiles into a single image
    combined_img = combine_tiles(tiles, tile_size)

    # Convert the image to a NumPy array
    img_array = np.array(combined_img)
    return img_array

def Esri_relief_images(aExtent, zoom_level):
    minx, maxx, miny, maxy = aExtent
    x_min, y_min, x_max, y_max = extent_to_tile_indices(minx, miny, maxx, maxy, zoom_level)
    tile_size = 256
    tiles = []
    for y in range(y_min, y_max + 1):
        row = []
        for x in range(x_min, x_max + 1):
            tile = fetch_esri_relif_tile(zoom_level, x, y)
            row.append(tile)
        tiles.append(row)

    # Combine the tiles into a single image
    combined_img = combine_tiles(tiles, tile_size)
    # Convert the image to a NumPy array
    img_array = np.array(combined_img)
    return img_array

def Esri_hydro_images(aExtent, zoom_level):
    # Calculate the tile coordinates for the domain
    minx, maxx, miny,  maxy = aExtent

    x_min, y_min, x_max, y_max = extent_to_tile_indices(minx, miny, maxx, maxy, zoom_level)
    tile_size = 256
    tiles = []
    for y in range(y_min, y_max + 1):
        row = []
        for x in range(x_min, x_max + 1):
            tile = fetch_esri_hydro_tile(zoom_level, x, y)
            row.append(tile)
        tiles.append(row)

    # Combine the tiles into a single image
    combined_img = combine_tiles(tiles, tile_size)
    # Convert the image to a NumPy array
    img_array = np.array(combined_img)
    return img_array

def calculate_zoom_level(scale_denominator, pProjection, dpi=96, tile_width=256, tile_height=256):
    """
    Calculates the appropriate zoom level based on the scale denominator, CRS, and DPI.

    Args:
      scale_denominator: The scale denominator of the map.
      pProjection: The coordinate reference system (CRS) of the map in WKT format.
      dpi: The dots per inch (DPI) of the display (default: 96).
      tile_width: The width of the tiles (default: 512).
      tile_height: The height of the tiles (default: 512).

    Returns:
      The calculated zoom level.
    """

    # Convert the projection to spatial reference
    pSpatial_reference_target = osr.SpatialReference()
    pSpatial_reference_target.ImportFromWkt(pProjection)
    meters_per_unit = pSpatial_reference_target.GetLinearUnits()

    # Calculate pixel span considering DPI
    pixel_size_in_meters = 0.00028
    pixel_span = scale_denominator * pixel_size_in_meters / meters_per_unit / (dpi / 96.0)

    # Calculate tile spans
    tile_span_x = tile_width * pixel_span
    tile_span_y = tile_height * pixel_span

    # Calculate zoom level based on tile spans
    zoom_level = int(math.log2(40075016.68557849 / max(tile_span_x, tile_span_y)))

    return zoom_level

def calculate_scale_denominator(domain_boundary, image_size, dpi=96):

    """
    Calculates the scale denominator for a map based on the domain boundary and desired image size.

    Args:
      domain_boundary: A list or tuple containing the minimum and maximum x and y coordinates of the domain.
      image_size: A tuple containing the desired width and height of the image in pixels.
      crs: The coordinate reference system (CRS) of the domain boundary (default: WGS84).

    Returns:
      The calculated scale denominator.
    """

    # Extract domain boundaries
    min_x, max_x, min_y, max_y = domain_boundary

    # Calculate domain width and height in map units
    domain_width = max_x - min_x
    domain_height = max_y - min_y

    # Get meters per unit for the CRS

    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    dLatitude_mean = (min_y + max_y) / 2.0
    meters_per_dgree = degree_to_meter(1.0, dLatitude_mean) #meter per degree
    meters_to_inches = 39.3701
    image_width_in_inches = image_size[0] / dpi
    image_height_in_inches = image_size[1] / dpi

    # Calculate scale denominator, meter per inch
    scale_denominator = np.max([domain_width, domain_height]) * meters_per_dgree * meters_to_inches \
            / np.max([image_width_in_inches, image_height_in_inches])

    return scale_denominator
