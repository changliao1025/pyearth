sFormat_y = '%.1E'

height = 1893434.0
print( '{:.1E}'.format(height)  )

print( sFormat_y % height )

format_e = '.1E'
format_s = '.2f'
txt = '{height:'+ format_e+'}​​​​​​​​​​​​​'
print(type(sFormat_y))

print(type(format_s))

print(type(txt))
print(txt.format(height = height))
import string


class MyFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        if format_spec == 'm':
            d= super().format_field(value, 'e').replace('e+', 'E')
            return d
        else:
            c =super().format_field(value, format_spec)
            return c


fmt = MyFormatter()
v = 1e+06
print(fmt.format('{:e}, {:m}', v, v))  # -> 1.000000e+06, 1.000000e06






