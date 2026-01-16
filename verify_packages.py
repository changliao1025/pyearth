from setuptools import find_packages

packages = find_packages()
print(f'Total packages found: {len(packages)}')

newly_added = [
    'pyearth.toolbox.analysis.difference',
    'pyearth.toolbox.analysis.intersect',
    'pyearth.toolbox.data.gsim',
    'pyearth.toolbox.data.list',
    'pyearth.toolbox.data.ocean',
    'pyearth.toolbox.data.vector',
    'pyearth.toolbox.mesh',
    'pyearth.toolbox.mesh.latlon',
    'pyearth.toolbox.mesh.hexaon',
    'pyearth.toolbox.mesh.square',
    'pyearth.toolbox.mesh.algorithm',
    'pyearth.toolbox.cython',
    'pyearth.visual.geovista',
    'pyearth.visual.ridgeplot'
]

print('\nNewly added packages verification:')
for p in newly_added:
    status = '✓' if p in packages else '✗'
    print(f'  {status} {p}')

all_present = all(p in packages for p in newly_added)
print(f'\nAll packages present: {all_present}')
