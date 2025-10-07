def clean_geometry(geometry, tolerance=1e-6, verbose=False):
    """
    Clean geometry artifacts and issues that can occur after union operations

    :param geometry: OGR geometry object to clean
    :param tolerance: tolerance for simplification and buffer operations
    :param verbose: whether to print cleaning information
    :return: cleaned geometry or None if cleaning failed
    """
    if geometry is None:
        return None

    try:
        original_type = geometry.GetGeometryName()

        # Step 1: Remove duplicate vertices and fix self-intersections
        if hasattr(geometry, 'MakeValid'):
            # Use MakeValid if available (GDAL 3.0+)
            geometry = geometry.MakeValid()
        else:
            # Fallback: use buffer trick to fix self-intersections
            geometry = geometry.Buffer(0)

        # Step 2: Simplify to remove collinear points and tiny segments
        if tolerance > 0:
            geometry = geometry.Simplify(tolerance)

        # Step 3: Remove slivers using small negative then positive buffer
        # This removes very thin areas that might appear as lines
        buffer_tolerance = tolerance * 10 if tolerance > 0 else 1e-8
        geometry = geometry.Buffer(-buffer_tolerance).Buffer(buffer_tolerance)

        # Step 4: Ensure we still have the correct geometry type
        if geometry and geometry.GetGeometryName() != original_type:
            if verbose:
                print(f"Warning: Geometry type changed from {original_type} to {geometry.GetGeometryName()} during cleaning")

        # Step 5: Final validation
        if geometry and not geometry.IsValid():
            if verbose:
                print("Warning: Geometry still invalid after cleaning, applying final buffer fix")
            geometry = geometry.Buffer(0)

        return geometry

    except Exception as e:
        if verbose:
            print(f"Warning: Failed to clean geometry: {str(e)}")
        return geometry  # Return original if cleaning fails