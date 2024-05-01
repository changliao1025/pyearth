import itertools
from matplotlib.patheffects import Stroke, Normal
import numpy as np
import cartopy.mpl.geoaxes

def zebra_frame(self, lw=3, crs=None, zorder=None, iFlag_outer_frame_in = None):    
    # Alternate black and white line segments
    bws = itertools.cycle(["k", "w"])
    self.spines["geo"].set_visible(False)
    
    if iFlag_outer_frame_in is not None:
        #get the map spatial reference        
        left, right, bottom, top = self.get_extent()
        crs_map = self.projection
        xticks = np.arange(left, right+(right-left)/9, (right-left)/8)
        yticks = np.arange(bottom, top+(top-bottom)/9, (top-bottom)/8)
        #check spatial reference are the same           
        pass
    else:        
        crs_map =  crs
        xticks = sorted([*self.get_xticks()])
        xticks = np.unique(np.array(xticks))        
        yticks = sorted([*self.get_yticks()])
        yticks = np.unique(np.array(yticks))        

    for ticks, which in zip([xticks, yticks], ["lon", "lat"]):
        for idx, (start, end) in enumerate(zip(ticks, ticks[1:])):
            bw = next(bws)
            if which == "lon":
                xs = [[start, end], [start, end]]
                ys = [[yticks[0], yticks[0]], [yticks[-1], yticks[-1]]]
            else:
                xs = [[xticks[0], xticks[0]], [xticks[-1], xticks[-1]]]
                ys = [[start, end], [start, end]]

            # For first and last lines, used the "projecting" effect
            capstyle = "butt" if idx not in (0, len(ticks) - 2) else "projecting"
            for (xx, yy) in zip(xs, ys):
                self.plot(xx, yy, color=bw, linewidth=max(0, lw - self.spines["geo"].get_linewidth()*2), clip_on=False,
                    transform=crs_map, zorder=zorder, solid_capstyle=capstyle,
                    # Add a black border to accentuate white segments
                    path_effects=[
                        Stroke(linewidth=lw, foreground="black"),
                        Normal(),
                    ],
                )

setattr(cartopy.mpl.geoaxes.GeoAxes, 'zebra_frame', zebra_frame)
