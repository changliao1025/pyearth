import cartopy.crs as ccrs
from osgeo import osr
def cartopy_define_projection(pProjection_in):
    
    km_proj = {'lon_0': 'central_longitude',
               'lat_0': 'central_latitude',
               'x_0': 'false_easting',
               'y_0': 'false_northing',
               'lat_ts': 'latitude_true_scale',
               'o_lon_p': 'central_rotated_longitude',
               'o_lat_p': 'pole_latitude',
               'sKey': 'scale_factor',
               'zone': 'zone',
               }
    km_globe = {'a': 'semimajor_axis',
                'b': 'semiminor_axis',
                }
    km_std = {'lat_1': 'lat_1',
              'lat_2': 'lat_2',
              }
    kw_proj = dict()
    kw_globe = dict()
    kw_std = dict()
    srs = pProjection_in
    
    if srs.IsProjected:
        sProj = srs.GetAttrValue('projcs')
        #get parameter
        print(sProj)
        
        srs_new = srs.ExportToProj4()
        for s in srs_new.split('+'):
            s = s.split('=')
            if len(s) != 2:
                continue
            sKey = s[0].strip()
            sValue = s[1].strip()
            print(sKey)
            try:
                sValue = float(sValue)
            except:
                pass
            if sKey == 'proj':
                if sValue == 'tmerc':
                    cl = ccrs.TransverseMercator
                if sValue == 'lcc':
                    cl = ccrs.LambertConformal
                if sValue == 'merc':
                    cl = ccrs.Mercator
                if sValue == 'utm':
                    cl = ccrs.UTM
                if sValue == 'stere':
                    cl = ccrs.Stereographic
                if sValue == 'ob_tran':
                    cl = ccrs.RotatedPole
                if sValue =='aea':
                    cl = ccrs.AlbersEqualArea

            if sKey in km_proj:
                kw_proj[km_proj[sKey]] = sValue
            if sKey in km_globe:
                kw_globe[km_globe[sKey]] = sValue
            if sKey in km_std:
                kw_std[km_std[sKey]] = sValue

        globe = None
        if kw_globe:
            globe = ccrs.Globe(ellipse='sphere', **kw_globe)
        if kw_std:
            kw_proj['standard_parallels'] = (kw_std['lat_1'], kw_std['lat_2'])

        # mercatoooor
        if cl.__name__ == 'Mercator':
            kw_proj.pop('false_easting', None)
            kw_proj.pop('false_northing', None)
            if LooseVersion(cartopy.__version__) < LooseVersion('0.15'):
                kw_proj.pop('latitude_true_scale', None)
        elif cl.__name__ == 'Stereographic':
            kw_proj.pop('scale_factor', None)
            if 'latitude_true_scale' in kw_proj:
                kw_proj['true_scale_latitude'] = kw_proj['latitude_true_scale']
                kw_proj.pop('latitude_true_scale', None)
        elif cl.__name__ == 'RotatedPole':
            if 'central_longitude' in kw_proj:
                kw_proj['pole_longitude'] = kw_proj['central_longitude'] - 180
                kw_proj.pop('central_longitude', None)
        elif cl.__name__ == 'AlbersEqualArea':
            
                
            #kw_proj.pop('central_longitude', None)
            #kw_proj.pop('false_easting', None)
            #kw_proj.pop('false_northing', None)
            pass
        else:
            kw_proj.pop('latitude_true_scale', None)

        return cl(globe=globe, **kw_proj)


    else:

        if srs.IsGeographic: 
            sGcs = srs.GetAttrValue('geogcs')
            print(sGcs)
        else:
                pass
             
    return sProjection