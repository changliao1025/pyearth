def log_formatter(x):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

def float_formatter(x):
    a = '{:.1f}'.format(x)
    return a