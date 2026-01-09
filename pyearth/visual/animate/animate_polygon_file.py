import os
import json
import numpy as np
from osgeo import osr, gdal, ogr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib import animation
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

pProjection = ccrs.PlateCarree()  # for latlon data only
iFigwidth_default = 8
iFigheight_default = 6


def animate_vector_polygon_data(
    sFilename_mesh_in,
    sFilename_animation_json_in,
    sFilename_animation_out,
    iFlag_type_in=None,
    iFigwidth_in=None,
    iFigheight_in=None,
    aExtent_in=None,
    sTitle_in=None,
    pProjection_map_in=None,
):

    if iFigwidth_in is None:
        iFigwidth_in = iFigwidth_default

    if iFigheight_in is None:
        iFigheight_in = iFigheight_default

    if iFlag_type_in is None:
        iFlag_type_in = 1

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180

    aCell_animation = list()
    aData = list()
    aData_raw = list()
    aLongitude = list()
    aLatitude = list()
    aAVertex = list()
    # get domain range
    pDriver = ogr.GetDriverByName("GeoJSON")
    pDataset = pDriver.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)  # WGS84 lat/lon
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type == "POLYGON":
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            dLon_max = np.max([dLon_max, np.max(aCoords_gcs[:, 0])])
            dLon_min = np.min([dLon_min, np.min(aCoords_gcs[:, 0])])
            dLat_max = np.max([dLat_max, np.max(aCoords_gcs[:, 1])])
            dLat_min = np.min([dLat_min, np.min(aCoords_gcs[:, 1])])

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(
            central_longitude=0.50 * (dLon_max + dLon_min),
            central_latitude=0.50 * (dLat_max + dLat_min),
            globe=None,
        )

    fig = plt.figure(dpi=100)
    fig.set_figwidth(iFigwidth_in)
    fig.set_figheight(iFigheight_in)
    ax = fig.add_axes([0.1, 0.1, 0.65, 0.8], projection=pProjection_map)
    ax.set_global()
    marginx = (dLon_max - dLon_min) / 20
    marginy = (dLat_max - dLat_min) / 20
    if aExtent_in is None:
        aExtent = [
            dLon_min - marginx,
            dLon_max + marginx,
            dLat_min - marginy,
            dLat_max + marginy,
        ]
    else:
        aExtent = aExtent_in

    ax.set_extent(aExtent)
    ax.coastlines()  # resolution='110m')
    if iFlag_type_in == 1:  # full
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=1,
            color="gray",
            alpha=0.3,
            linestyle="--",
        )
        gl.xlabel_style = {"size": 8, "color": "k", "rotation": 0, "ha": "right"}
        gl.ylabel_style = {"size": 8, "color": "k", "rotation": 90, "weight": "normal"}
    else:  # track mode
        pass

    cmap = cm.get_cmap("Spectral")
    if sTitle_in is not None:
        # setting a title for the plot
        # sText = "Priority flood in HexWatershed"
        ax.text(
            0.5,
            1.08,
            sTitle_in,
            verticalalignment="center",
            horizontalalignment="center",
            transform=ax.transAxes,
            color="black",
            fontsize=9,
        )
    cmap_reversed = cmap.reversed()

    # get dataset
    with open(sFilename_animation_json_in) as json_file:
        aCell_animation = json.load(json_file)
        ncell_animation = len(aCell_animation)

    for i in range(ncell_animation):
        pCell_animation = aCell_animation[i]
        aData.append(float(pCell_animation["dElevation"]))
        aData_raw.append(float(pCell_animation["dElevation_raw"]))
        aLongitude.append(float(pCell_animation["dLongitude_center_degree"]))
        aLatitude.append(float(pCell_animation["dLatitude_center_degree"]))
        aAVertex.append(pCell_animation["vVertex"])

    # convert to numpy array
    aData = np.array(aData)
    aData_raw = np.array(aData_raw)
    aLongitude = np.array(aLongitude)
    aLatitude = np.array(aLatitude)
    dData_max = np.max(aData)
    dData_min = np.min(aData)

    sText = ""
    x1 = 0.0
    y1 = 0.0
    pArtist1 = ax.text(
        x1,
        y1,
        sText,
        verticalalignment="center",
        horizontalalignment="right",
        transform=ax.transAxes,
        color="black",
        fontsize=12,
    )
    x2 = 0.0
    y2 = 0.0
    pArtist2 = ax.text(
        x2,
        y2,
        sText,
        verticalalignment="center",
        horizontalalignment="left",
        transform=ax.transAxes,
        color="black",
        fontsize=12,
    )

    # trasform data
    norm = plt.Normalize(dData_min, dData_max)
    sm = plt.cm.ScalarMappable(cmap=cmap_reversed, norm=norm)
    sm.set_array(aData)
    fig.canvas.draw()
    # Section 2
    ax_pos = ax.get_position()  # get the original position
    # use this ax to set the colorbar ax position
    ax_cb = fig.add_axes([ax_pos.x1 + 0.06, ax_pos.y0, 0.02, ax_pos.height])
    cb = fig.colorbar(sm, cax=ax_cb)
    sUnit = r"Unit: m"
    cb.ax.get_yaxis().set_ticks_position("right")
    cb.ax.get_yaxis().labelpad = 5
    cb.ax.set_ylabel(sUnit, rotation=90)
    cb.ax.get_yaxis().set_label_position("left")
    cb.ax.tick_params(labelsize=6)

    # calculate ahead
    aColor_index = (aData - dData_min) / (dData_max - dData_min)
    aColors = [cmap_reversed(aColor_index[i]) for i in range(ncell_animation)]

    # Precompute polygon locations and add them to the axes
    aPolygon = [
        mpatches.Polygon(
            np.array(
                [
                    [vertex["dLongitude_degree"], vertex["dLatitude_degree"]]
                    for vertex in aVertex
                ]
            ),
            closed=True,
            transform=ccrs.Geodetic(),
            facecolor="none",  # Initially invisible
            edgecolor="none",  # Optional: no border
        )
        for aVertex in aAVertex
    ]

    # Add polygons to the axes
    for polygon in aPolygon:
        ax.add_patch(polygon)

    # Initialize artists
    pArtist0 = ax.add_patch(aPolygon[0])
    pArtist1 = ax.text(
        0.0,
        0.0,
        "",
        verticalalignment="center",
        horizontalalignment="right",
        transform=ax.transAxes,
        fontsize=12,
    )
    pArtist2 = ax.text(
        0.0,
        0.0,
        "",
        verticalalignment="center",
        horizontalalignment="left",
        transform=ax.transAxes,
        fontsize=12,
    )

    def animate(i):
        rgb = aColors[i]
        aPolygon[i].set_facecolor(rgb)
        pArtist0.set_xy(aPolygon[i].get_xy())
        x = (aLongitude[i] - dLon_min) / (dLon_max - dLon_min)
        y = (aLatitude[i] - dLat_min) / (dLat_max - dLat_min)
        dummy0, dummy = aData_raw[i], aData[i]
        y1, y2 = (y + 0.05, y - 0.05) if dummy0 > dummy else (y - 0.05, y + 0.05)
        pArtist1.set_position((x, y1))
        pArtist1.set_text(f"Before: {dummy0:.2f}m")
        pArtist2.set_position((x, y2))
        pArtist2.set_text(f"After: {dummy:.2f}m")

        return pArtist0, pArtist1, pArtist2

    # remove if already exist
    # check output file extension
    if os.path.exists(sFilename_animation_out):
        os.remove(sFilename_animation_out)
    if sFilename_animation_out.endswith(".gif"):
        plt.rcParams["animation.convert_path"] = (
            "/share/apps/ImageMagick/7.1.0-52/bin/convert"
        )
        writer = "imagemagick"
        anim = animation.FuncAnimation(
            fig, animate, frames=ncell_animation, interval=400, blit=False
        )
        anim.save(sFilename_animation_out, writer=writer)

    else:  # use video instead
        plt.rcParams["animation.ffmpeg_path"] = (
            "/people/liao313/.conda/envs/pyflowline/bin/ffmpeg"  #'/share/apps/ffmpeg/6.1.1/bin/ffmpeg
        )
        # Writer = animation.writers['ffmpeg']
        # writer = Writer(fps=5, metadata=dict(artist='Chang Liao'), bitrate=-1)
        writer = animation.FFMpegWriter(
            fps=5, metadata=dict(artist="Chang Liao"), bitrate=-1
        )
        anim = animation.FuncAnimation(
            fig, animate, frames=ncell_animation, interval=1, blit=True
        )
        anim.save(sFilename_animation_out, writer=writer)
    return
