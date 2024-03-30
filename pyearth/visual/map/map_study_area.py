# plot a map of study area which has google earth image, slope, and dem
import numpy as np
from osgeo import osr, gdal, ogr

import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy as cpl

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.gis.spatialref.reproject_coodinates import reproject_coordinates, reproject_coordinates_batch
from pyearth.gis.location.Google_MetersPerPixel import Google_MetersPerPixel
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates


def plot_study_area(sFilename_dem_in,
                    sFilename_out,
                    sFilename_polygon_in=None,
                    sFilename_slope_in=None,
                    sFilename_point_in=None,
                    sFilename_polyline_in=None):

    try:
        import requests
    except ImportError as e:
        raise ImportError(
            "The package 'requests' is required for this function to run.") from e

    # set up figure size and dpi
    fig = plt.figure(dpi=300)
    fig.set_figwidth(4)
    fig.set_figheight(4)
    # ==============================================================================
    # plot dem first
    # ==============================================================================

    dummy = gdal_read_geotiff_file(sFilename_dem_in)
    aImage_in = dummy['dataOut']
    dResolution_x = dummy['pixelWidth']
    dResolution_y = dummy['pixelHeight']
    dOriginX = dummy['originX']
    dOriginY = dummy['originY']
    nrow = dummy['nrow']
    ncolumn = dummy['ncolumn']
    missing_value = dummy['missingValue']
    pProjection = dummy['projection']
    pSpatialRef_source = dummy['spatialReference']

    # set up dem projection
    dem_proj = cpl.crs.AlbersEqualAre(central_longitude=pSpatialRef_source.GetProjPar('longitude_of_center'),
                                   central_latitude=pSpatialRef_source.GetProjPar(
                                       'latitude_of_center'),
                                   standard_parallels=(pSpatialRef_source.GetProjPar('standard_parallel_1'),
                                                       pSpatialRef_source.GetProjPar('standard_parallel_2')))

    # set up dem plot
    ax_dem = fig.add_axes([0.3, 0.15, 0.6, 0.55], projection=dem_proj)
    ax_dem.set_xmargin(0.05)
    ax_dem.set_ymargin(0.10)

    dLon_min = dOriginX
    dLon_max = dOriginX + ncolumn * dResolution_x
    dLat_max = dOriginY
    dLat_min = dOriginY + nrow * dResolution_y

    pSpatialRef_target = osr.SpatialReference()
    pSpatialRef_target.ImportFromEPSG(4326)
    aLon = list()
    aLat = list()
    aLon.append(dLon_min)
    aLon.append(dLon_max)
    aLat.append(dLat_min)
    aLat.append(dLat_max)
    aLon, aLat = reproject_coordinates_batch(
        aLon, aLat, pSpatialRef_source, pSpatialRef_target)
    dLongitude_center = np.mean(aLon)
    dLatitude_center = np.mean(aLat)
    aImage_extent = [dLon_min - dResolution_x, 
                     dLon_max + dResolution_x, 
                     dLat_min + dResolution_y,  
                     dLat_max - dResolution_y]

    aImage_in[np.where(aImage_in == missing_value)] = np.nan
    #
    demplot = ax_dem.imshow(aImage_in, origin='upper',
                            extent=aImage_extent, cmap=mpl.cm.terrain, transform=dem_proj)
    # change all spines
    for axis in ['top', 'bottom', 'left', 'right']:
        ax_dem.spines[axis].set_linewidth(0.5)
        # increase tick width
    ax_dem.tick_params(width=0.5)

    # plot flowline if available
    pDriver = ogr.GetDriverByName('GeoJSON')
    if sFilename_polyline_in is None:
        pass
    else:
        pDataset = pDriver.Open(sFilename_polyline_in, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type == 'LINESTRING':
                # dummy0 = loads( pGeometry_in.ExportToWkt() )
                # aCoords_gcs = dummy0.coords
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = np.array(aCoords_gcs)
                aCoords_gcs = aCoords_gcs[:, 0:2]
                nvertex = len(aCoords_gcs)

                codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int)
                codes[0] = mpl.path.Path.MOVETO
                path = mpl.path.Path(aCoords_gcs, codes)
                x, y = zip(*path.vertices)
                line, = ax_dem.plot(
                    x, y, color='b', transform=cpl.crs.PlateCarree())
                lID = lID + 1

    # plot point if available
    if sFilename_point_in is None:
        pass
    else:
        pDataset = pDriver.Open(sFilename_point_in, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type == 'POINT':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = np.array(aCoords_gcs)
                aCoords_gcs = aCoords_gcs[:, 0:2]
                nvertex = len(aCoords_gcs)

                codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int)
                codes[0] = mpl.path.Path.MOVETO
                path = mpl.path.Path(aCoords_gcs, codes)
                x, y = zip(*path.vertices)
                line, = ax_dem.plot(
                    x, y, color='b', transform=cpl.crs.PlateCarree())
                lID = lID + 1

    # draw gridlines
    gl = ax_dem.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.3, linestyle='--')

    gl.left_labels = False
    gl.top_labels = False
    XTEXT_SIZE = 8
    YTEXT_SIZE = 8
    # to facilitate text rotation at bottom edge, ...
    # text justification: 'ha':'right' is used to avoid clashing with map'sboundary
    # default of 'ha' is center, often causes trouble when text rotation isnot zero
    gl.xlabel_style = {'size': XTEXT_SIZE,
                       'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': YTEXT_SIZE, 'color': 'k',
                       'rotation': 90, 'weight': 'normal'}

    # draw colorbar
    ax_cb = fig.add_axes([0.2, 0.2, 0.02, 0.5])
    cb = plt.colorbar(demplot, cax=ax_cb, extend='both')
    cb.ax.get_yaxis().set_ticks_position('left')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel('Unit: meter', rotation=270)
    cb.ax.tick_params(labelsize=6)

    # ==============================================
    # google earth
    ge_proj = cpl.crs.Mercator()  # central_longitude=dLongitude_center)

    pSpatialRef_source = osr.SpatialReference()
    pSpatialRef_source.ImportFromEPSG(4326)
    pSpatialRef_target = osr.SpatialReference()

    iZoom = 7
    pSpatialRef_target.ImportFromWkt(ge_proj.to_wkt())
    uv_xcenter, uv_ycenter = reproject_coordinates(
        dLongitude_center, dLatitude_center, pSpatialRef_source, pSpatialRef_target)
    xsize_ge = 600
    ysize_ge = 600
    scale = Google_MetersPerPixel(iZoom)
    xrange = (uv_xcenter - (xsize_ge/2.0*scale),
              uv_xcenter + (xsize_ge/20*scale))
    yrange = (uv_ycenter - (ysize_ge/2.0*scale),
              uv_ycenter + (ysize_ge/20*scale))
    aImage_extent = [xrange[0], xrange[1], yrange[0], yrange[1]]
    ax_ge = fig.add_axes([0.1, 0.75, 0.2, 0.2], projection=ge_proj)

    sLongitude_center = "{:0f}".format(dLongitude_center)
    sLatitude_center = "{:0f}".format(dLatitude_center)
    sZoom = "{:0d}".format(iZoom)

    sResolution = "{:0d}".format(xsize_ge) + 'x' + "{:0d}".format(ysize_ge)
    sMap_type = "hybrid"
    #obtain GCP APIs from your system
    sGoogleMapAPI = ''
    sGoogleMap = "http://maps.googleapis.com/maps/api/staticmap?" + \
        "center=" + sLatitude_center + ',' + sLongitude_center + \
        "&zoom=" + sZoom + "&size=" + sResolution + \
        "&maptype="+sMap_type+"&sensor=false&format=png32" + \
        "key=" + sGoogleMapAPI
    r = requests.get(sGoogleMap)
    # wb mode is stand for write binary mode
    f = open('google_map.png', 'wb')
    # r.content gives content,
    # in this case gives image
    f.write(r.content)
    # close method of file object
    # save and close the file
    f.close()

    img = mpl.image.imread('google_map.png')
    ax_ge.imshow(img, extent=aImage_extent,       transform=ge_proj)
    gl_ge = ax_ge.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                            linewidth=1, color='gray', alpha=0.3, linestyle='--')
    # gl_ge.xlocator = mpl.ticker.FixedLocator(aLon_range)
    # gl_ge.ylocator = mpl.ticker.FixedLocator(aLat_range)
    gl_ge.xformatter = LongitudeFormatter()
    gl_ge.yformatter = LatitudeFormatter()
    gl_ge.right_labels = False
    gl_ge.bottom_labels = False

    # declare text size
    XTEXT_SIZE = 2.5
    YTEXT_SIZE = 2.5
    gl_ge.xlabel_style = {'size': XTEXT_SIZE,
                          'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl_ge.ylabel_style = {'size': YTEXT_SIZE,
                          'color': 'k', 'rotation': 90, 'weight': 'normal'}
    # add boundary
    if sFilename_polygon_in is None:
        pass
    else:
        pDataset = pDriver.Open(sFilename_polygon_in, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type == 'POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aCoords_gcs = np.array(aCoords_gcs)
                polygon = mpl.patches.Polygon(aCoords_gcs[:, 0:2], closed=True, linewidth=0.25,
                                           alpha=0.8, edgecolor='b',
                                           transform=cpl.crs.PlateCarree())
                ax_ge.add_patch(polygon)

    # histogram
    ax_histo = fig.add_axes([0.4, 0.8, 0.4, 0.15])
    if sFilename_slope_in is None:
        pass
    else:
        dummy = gdal_read_geotiff_file(sFilename_slope_in)
        aSlope = dummy[0]
        missing_value = dummy[7]
        aSlope = aSlope[np.where(aSlope != missing_value)]
        dMax_x = 60
        dMin_x = 0
        dSpace_x = 4
        ax_histo.hist(aSlope,  int((dMax_x-dMin_x) / dSpace_x),
                      color="skyblue", ec="skyblue")
        sLabel_x = 'Slope (degree)'
        sLabel_y = 'Frequency'
        ax_histo.set_xlabel(sLabel_x, fontsize=4)
        ax_histo.set_ylabel(sLabel_y, fontsize=4)
        ax_histo.set_xlim(dMin_x, dMax_x)
        formatter = mpl.ticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        ax_histo.yaxis.set_major_formatter(formatter)
        ax_histo.tick_params(axis='x', labelsize=5)
        ax_histo.tick_params(axis='y', labelsize=5)
        ax_histo.yaxis.offsetText.set_fontsize(5)

    plt.savefig(sFilename_out, bbox_inches='tight')

    return
