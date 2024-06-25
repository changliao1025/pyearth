def pick_colormap_hydrology(sVariable_in):
  """
  Selects a colormap appropriate for a given hydrologic sVariable_in.

  Args:
    sVariable_in: The name of the hydrologic sVariable_in (e.g., 'precipitation', 'streamflow', 'soil moisture').

  Returns:
    A string representing the name of a matplotlib colormap.
  """

  colormap_dict = {
      'precipitation': 'Blues',
      'streamflow': 'bwr',
      'soil moisture': 'RdGy',
      'evapotranspiration': 'Greens',
      'water depth': 'Blues',
      'snow water equivalent': 'Greys',
      'channel width': 'YlGnBu',
      'channel depth': 'Blues',
      'flood fraction': 'Blues',
  }

  return colormap_dict.get(sVariable_in.lower(), 'viridis')

