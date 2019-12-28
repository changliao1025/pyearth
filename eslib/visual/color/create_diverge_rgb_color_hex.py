def create_diverge_rgb_color_hex(ncolor):

    if ncolor < 3 or ncolor > 12:
        return -1
    else:
        if ncolor == 3:
            pass
        else:
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
                pass
        return colors_hex