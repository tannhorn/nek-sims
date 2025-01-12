import gmsh


def generate_points(
    coord_x: list[float], coord_y: list[float]
) -> tuple[list[int], dict[tuple[int, int], int]]:
    """
    Generates points at each combination of coordinates from the input lists in Gmsh.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.

    Returns:
    tuple: A tuple containing a list of Gmsh point tags and a mapping of coordinates to point tags.
    """
    # Keep track of the generated point tags
    point_tags = []
    point_map = {}

    # Generate points at each combination of coordinates
    for i, x in enumerate(coord_x):
        for j, y in enumerate(coord_y):
            tag = gmsh.model.geo.addPoint(x, y, 0)  # Add a point at (x, y, 0)
            point_tags.append(tag)
            point_map[(i, j)] = tag

    return point_tags, point_map


def generate_lines(
    coord_x: list[float], coord_y: list[float], point_map: dict[tuple[int, int], int]
) -> dict[tuple[int, int, int, int], int]:
    """
    Generates lines between adjacent points based on the given point map.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    point_map (dict[tuple[int, int], int]): Mapping of coordinates to point tags.

    Returns:
    dict[tuple[int, int, int, int], int]: A dictionary mapping point indices to Gmsh line tags.
    """
    line_map = {}

    # Generate lines parallel to x-direction
    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y)):
            start_tag = point_map[(i, j)]
            end_tag = point_map[(i + 1, j)]
            line_map[(i, j, i + 1, j)] = gmsh.model.geo.addLine(start_tag, end_tag)

    # Generate lines parallel to y-direction
    for i in range(len(coord_x)):
        for j in range(len(coord_y) - 1):
            start_tag = point_map[(i, j)]
            end_tag = point_map[(i, j + 1)]
            line_map[(i, j, i, j + 1)] = gmsh.model.geo.addLine(start_tag, end_tag)

    return line_map


def apply_transfinite_curves(
    lines: dict[tuple[int, int, int, int], int],
    coord_x: list[float],
    coord_y: list[float],
    disc_x: list[int],
    disc_y: list[int],
    prog_x: list[float],
    prog_y: list[float],
):
    """
    Applies transfinite curves to the generated lines based on the provided progressions.

    Parameters:
    lines (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    disc_x (list[int]): Discretization values for lines parallel to the x-direction.
    disc_y (list[int]): Discretization values for lines parallel to the y-direction.
    prog_x (list[float]): Progression values for lines parallel to the x-direction.
    prog_y (list[float]): Progression values for lines parallel to the y-direction.
    """
    # Apply transfinite curves to lines parallel to x-direction
    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y)):
            line_tag = lines[(i, j, i + 1, j)]
            gmsh.model.geo.mesh.setTransfiniteCurve(
                line_tag, disc_x[i], "Progression", prog_x[i]
            )

    # Apply transfinite curves to lines parallel to y-direction
    for i in range(len(coord_x)):
        for j in range(len(coord_y) - 1):
            line_tag = lines[(i, j, i, j + 1)]
            gmsh.model.geo.mesh.setTransfiniteCurve(
                line_tag, disc_y[j], "Progression", prog_y[j]
            )


def generate_surfaces(
    coord_x: list[float],
    coord_y: list[float],
    point_map: dict[tuple[int, int], int],
    line_map: dict[tuple[int, int, int, int], int],
) -> list[int]:
    """
    Generates transfinite surfaces for each rectangle enclosed by 4 lines.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    point_map (dict[tuple[int, int], int]): Mapping of coordinates to point tags.
    line_map (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.

    Returns:
    list[int]: A list of Gmsh surface tags.
    """
    surface_tags = []

    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y) - 1):
            # Retrieve the 4 lines that define the boundary of the rectangle
            l1 = line_map[(i, j, i + 1, j)]
            l2 = line_map[(i + 1, j, i + 1, j + 1)]

            l3 = line_map[(i, j + 1, i + 1, j + 1)]
            l4 = line_map[(i, j, i, j + 1)]

            # Create a curve loop and surface
            curve_loop = gmsh.model.geo.addCurveLoop(
                [l1, l2, -l3, -l4]
            )  # Negative sense for l3 and l4
            surface = gmsh.model.geo.addPlaneSurface([curve_loop])

            # Make the surface transfinite
            gmsh.model.geo.mesh.setTransfiniteSurface(surface)
            surface_tags.append(surface)

    return surface_tags


def add_physical_groups(
    coord_x: list[float],
    coord_y: list[float],
    lines: dict[tuple[int, int, int, int], int],
    surfaces: list[int],
):
    """
    Adds distinct physical groups to boundary lines and the entire surface.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    lines (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.
    surfaces (list[int]): List of Gmsh surface tags.
    """
    # Physical group for lines at xmin (left boundary)
    left_lines = [lines[(0, j, 0, j + 1)] for j in range(len(coord_y) - 1)]
    gmsh.model.addPhysicalGroup(1, left_lines, tag=1)

    # Physical group for lines at xmax (right boundary)
    right_lines = [
        lines[(len(coord_x) - 1, j, len(coord_x) - 1, j + 1)]
        for j in range(len(coord_y) - 1)
    ]
    gmsh.model.addPhysicalGroup(1, right_lines, tag=2)

    # Physical group for lines at ymin (bottom boundary)
    bottom_lines = [lines[(i, 0, i + 1, 0)] for i in range(len(coord_x) - 1)]
    gmsh.model.addPhysicalGroup(1, bottom_lines, tag=3)

    # Physical group for lines at ymax (top boundary)
    top_lines = [
        lines[(i, len(coord_y) - 1, i + 1, len(coord_y) - 1)]
        for i in range(len(coord_x) - 1)
    ]
    gmsh.model.addPhysicalGroup(1, top_lines, tag=4)

    # Physical group for all surfaces
    gmsh.model.addPhysicalGroup(2, surfaces, tag=5)


def mesh_and_save(surfaces: list[int], filename: str):
    """
    Meshes the surfaces, recombines the mesh into quads, and writes to a .msh file.

    Parameters:
    surfaces (list[int]): List of Gmsh surface tags.
    filename (str): Output .msh file name.
    """
    # Mesh the surfaces
    for surface in surfaces:
        gmsh.model.geo.mesh.setRecombine(2, surface)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Save the mesh to a file
    gmsh.write(filename)


# Example usage
if __name__ == "__main__":
    gmsh.initialize()
    gmsh.model.add("2D Mesh Points and Lines")

    coord_x = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
    coord_y = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.0]
    disc_x = [2, 3, 4, 5, 6]  # Discretization for x-direction lines
    disc_y = [2, 3, 4, 5, 6]  # Discretization for y-direction lines
    prog_x = [0.9, 1.1, 1.0, 0.7, 1.2]  # Progression for x-direction lines
    prog_y = [1, 1, 1, 1, 1]  # Progression for y-direction lines

    points, point_map = generate_points(coord_x, coord_y)
    print(f"Generated {len(points)} points: {points}")

    line_map = generate_lines(coord_x, coord_y, point_map)
    print(f"Generated {len(line_map)} lines.")

    apply_transfinite_curves(line_map, coord_x, coord_y, disc_x, disc_y, prog_x, prog_y)
    print(f"Meshed {len(line_map)} lines.")

    surfaces = generate_surfaces(coord_x, coord_y, point_map, line_map)
    print(f"Generated {len(surfaces)} surfaces: {surfaces}")

    add_physical_groups(coord_x, coord_y, line_map, surfaces)
    print(f"Added physical groups.")

    mesh_and_save(surfaces, "cross.msh")

    # Finalize Gmsh
    gmsh.finalize()
