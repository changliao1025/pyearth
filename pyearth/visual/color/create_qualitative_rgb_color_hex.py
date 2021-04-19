def create_qualitative_rgb_color_hex(ncolor, iFlag_reverse_in = None):
    """
    choose diverge color from: https://colorbrewer2.org/
    """
    if ncolor < 3 or ncolor > 12:
        return -1
    else:
        if ncolor == 3:
            colors_hex = [  '#fc8d59',\
                                '#ffffbf',\
                                '#91bfdb']

        else:
            if ncolor == 4:
                colors_hex = [  '#d7191c',\
                                '#fdae61',\
                                '#abdda4',\
                                '#2b83ba']
            if ncolor == 5:
                colors_hex = [  '#d7191c',\
                                '#fdae61',\
                                '#abdda4',\
                                '#2b83ba']
            if ncolor == 6:
                colors_hex = [  '#66c2a5',\
                                '#fc8d62',\
                                '#8da0cb',\
                                '#e78ac3',\
                                '#a6d854',\
                                '#ffd92f' ]


                colors_hex = [  '#e41a1c',\
                                '#377eb8',\
                                '#4daf4a',\
                                '#984ea3',\
                                '#ff7f00',\
                                '#ffff33']
            if ncolor == 7:
                colors_hex = [  '#e41a1c',\
                                '#377eb8',\
                                '#4daf4a',\
                                '#984ea3',\
                                '#ff7f00',\
                                '#ffff33',\
                                '#a65628']


                colors_hex = [  '#1b9e77',\
                                '#d95f02',\
                                '#7570b3',\
                                '#e7298a',\
                                '#66a61e',\
                                '#e6ab02',\
                                '#a6761d']
            if ncolor == 8:
                colors_hex = [  '#e41a1c' ,\
                                '#377eb8' ,\
                                '#4daf4a' ,\
                                '#984ea3' ,\
                                '#ff7f00' ,\
                                '#ffff33' ,\
                                '#a65628' ,\
                                '#f781bf' ]

            if ncolor == 9:
                colors_hex = [  '#a6cee3',\
                                '#1f78b4',\
                                '#b2df8a',\
                                '#33a02c',\
                                '#fb9a99',\
                                '#e31a1c',\
                                '#fdbf6f',\
                                '#ff7f00',\
                                '#cab2d6']
            if ncolor == 10:
                colors_hex = [  '#9e0142', \
                                '#d53e4f', \
                                '#f46d43', \
                                '#fdae61', \
                                '#fee08b', \
                                '#e6f598', \
                                '#abdda4', \
                                '#66c2a5', \
                                '#3288bd', \
                                '#5e4fa2']
            else:
                print('Too many colors are requested!')
                return -1
                
        #add the reverse feature
        if iFlag_reverse_in is not None:
            colors_hex.reverse()
        else:
            pass
        return colors_hex