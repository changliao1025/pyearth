import os
import cartopy.io.img_tiles as cimgt
class StadiaStamen(cimgt.Stamen):
    def _image_url(self, tile):
        x,y,z = tile
        #get sStadia_API_KEY from the environment variable
        sStadia_API_KEY = os.environ.get('STADIA_API_KEY')
        url = f"https://tiles.stadiamaps.com/tiles/stamen_terrain/{z}/{x}/{y}.png?api_key={sStadia_API_KEY}"
        return url