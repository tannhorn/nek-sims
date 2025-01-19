"""Prepare the Nek5000 mesh."""

import gmsh
import numpy as np

ZVAL: float = 0
LAMINAR: bool = True
RE: float = 1000


def generate_points(
    coord_x: list[float],
    coord_y: list[float],
    exclude: list[tuple[int, int]],
) -> tuple[list[int], dict[tuple[int, int], int]]:
    """
    Generates points at each combination of coordinates from the input lists in Gmsh.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.

    Returns:
    tuple: A tuple containing a list of Gmsh point tags and a mapping of coordinates to point tags.
    """
    # Keep track of the generated point tags
    point_tags = []
    point_map = {}

    # Generate points at each combination of coordinates
    for i, x in enumerate(coord_x):
        for j, y in enumerate(coord_y):
            if (i, j) in exclude:
                continue
            tag = gmsh.model.geo.addPoint(x, y, ZVAL)
            point_tags.append(tag)
            point_map[(i, j)] = tag

    return point_tags, point_map


def generate_lines(
    coord_x: list[float],
    coord_y: list[float],
    point_map: dict[tuple[int, int], int],
    exclude: list[tuple[int, int]],
) -> dict[tuple[int, int, int, int], int]:
    """
    Generates lines between adjacent points based on the given point map.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    point_map (dict[tuple[int, int], int]): Mapping of coordinates to point tags.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.

    Returns:
    dict[tuple[int, int, int, int], int]: A dictionary mapping point indices to Gmsh line tags.
    """
    line_map = {}

    # Generate lines parallel to x-direction
    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y)):
            if (i, j) in exclude or (i + 1, j) in exclude:
                continue
            start_tag = point_map[(i, j)]
            end_tag = point_map[(i + 1, j)]
            line_map[(i, j, i + 1, j)] = gmsh.model.geo.addLine(start_tag, end_tag)

    # Generate lines parallel to y-direction
    for i in range(len(coord_x)):
        for j in range(len(coord_y) - 1):
            if (i, j) in exclude or (i, j + 1) in exclude:
                continue
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
    exclude: list[tuple[int, int]],
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
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.
    """
    # Apply transfinite curves to lines parallel to x-direction
    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y)):
            if (i, j) in exclude or (i + 1, j) in exclude:
                continue
            line_tag = lines[(i, j, i + 1, j)]
            gmsh.model.geo.mesh.setTransfiniteCurve(
                line_tag, disc_x[i], "Progression", prog_x[i]
            )

    # Apply transfinite curves to lines parallel to y-direction
    for i in range(len(coord_x)):
        for j in range(len(coord_y) - 1):
            if (i, j) in exclude or (i, j + 1) in exclude:
                continue
            line_tag = lines[(i, j, i, j + 1)]
            gmsh.model.geo.mesh.setTransfiniteCurve(
                line_tag, disc_y[j], "Progression", prog_y[j]
            )


def generate_surfaces(
    coord_x: list[float],
    coord_y: list[float],
    line_map: dict[tuple[int, int, int, int], int],
    exclude: list[tuple[int, int]],
) -> list[int]:
    """
    Generates transfinite surfaces for each rectangle enclosed by 4 lines.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    point_map (dict[tuple[int, int], int]): Mapping of coordinates to point tags.
    line_map (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.

    Returns:
    list[int]: A list of Gmsh surface tags.
    """
    surface_tags = []

    for i in range(len(coord_x) - 1):
        for j in range(len(coord_y) - 1):
            if (
                (i, j) in exclude
                or (i + 1, j) in exclude
                or (i, j + 1) in exclude
                or (i + 1, j + 1) in exclude
            ):
                continue

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


def add_special_boundary_group(
    lines: dict[tuple[int, int, int, int], int],
    exclude: list[tuple[int, int]],
    existing_groups: list[list[int]],
    tag: int,
):
    """
    Adds lines bordering excluded points to a special physical group.

    Parameters:
    lines (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.
    existing_groups (list[list[int]]): List of line tags in existing physical groups.
    tag (int): The tag for the new physical group.
    """
    # Flatten the list of existing groups into a set for fast lookups
    existing_lines = set(line for group in existing_groups for line in group)

    special_lines = []

    for line_key, line_tag in lines.items():
        # Skip lines already in an existing physical group
        if line_tag in existing_lines:
            continue

        i1, j1, i2, j2 = line_key

        # Check adjacency to excluded points
        if i1 == i2:  # y-parallel line
            if (
                (i1 - 1, j1) in exclude
                or (i1 + 1, j1) in exclude
                or (i1 - 1, j2) in exclude
                or (i1 + 1, j2) in exclude
            ):
                special_lines.append(line_tag)

        elif j1 == j2:  # x-parallel line
            if (
                (i1, j1 - 1) in exclude
                or (i1, j1 + 1) in exclude
                or (i2, j2 - 1) in exclude
                or (i2, j2 + 1) in exclude
            ):
                special_lines.append(line_tag)

    # Add the special lines as a new physical group
    if special_lines:
        special_group = gmsh.model.addPhysicalGroup(1, special_lines, tag=tag)
        gmsh.model.setPhysicalName(1, special_group, "InternalSide")


def add_physical_groups(
    coord_x: list[float],
    coord_y: list[float],
    lines: dict[tuple[int, int, int, int], int],
    surfaces: list[int],
    exclude: list[tuple[int, int]],
):
    """
    Adds distinct physical groups to boundary lines and the entire surface.

    Parameters:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    lines (dict[tuple[int, int, int, int], int]): Mapping of point indices to Gmsh line tags.
    surfaces (list[int]): List of Gmsh surface tags.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.
    """
    # Physical group for lines at xmin (left boundary)
    left_lines = [
        lines[(0, j, 0, j + 1)]
        for j in range(len(coord_y) - 1)
        if (0, j) not in exclude and (0, j + 1) not in exclude
    ]
    left_group = gmsh.model.addPhysicalGroup(1, left_lines, tag=1)
    gmsh.model.setPhysicalName(1, left_group, "LeftSide")

    # Physical group for lines at xmax (right boundary)
    right_lines = [
        lines[(len(coord_x) - 1, j, len(coord_x) - 1, j + 1)]
        for j in range(len(coord_y) - 1)
        if (len(coord_x) - 1, j) not in exclude
        and (len(coord_x) - 1, j + 1) not in exclude
    ]
    right_group = gmsh.model.addPhysicalGroup(1, right_lines, tag=2)
    gmsh.model.setPhysicalName(1, right_group, "RightSide")

    # Physical group for lines at ymin (bottom boundary)
    bottom_lines = [
        lines[(i, 0, i + 1, 0)]
        for i in range(len(coord_x) - 1)
        if (i, 0) not in exclude and (i + 1, 0) not in exclude
    ]
    bottom_group = gmsh.model.addPhysicalGroup(1, bottom_lines, tag=3)
    gmsh.model.setPhysicalName(1, bottom_group, "BottomSide")

    # Physical group for lines at ymax (top boundary)
    top_lines = [
        lines[(i, len(coord_y) - 1, i + 1, len(coord_y) - 1)]
        for i in range(len(coord_x) - 1)
        if (i, len(coord_y) - 1) not in exclude
        and (i + 1, len(coord_y) - 1) not in exclude
    ]
    top_group = gmsh.model.addPhysicalGroup(1, top_lines, tag=4)
    gmsh.model.setPhysicalName(1, top_group, "TopSide")

    # Inner surfaces
    add_special_boundary_group(
        lines, exclude, [left_lines, right_lines, bottom_lines, top_lines], tag=5
    )

    # Physical group for all surfaces
    surface_group = gmsh.model.addPhysicalGroup(2, surfaces, tag=1000)
    gmsh.model.setPhysicalName(2, surface_group, "Domain")


def mesh_and_save(surfaces: list[int], filename: str):
    """
    Meshes the surfaces, recombines the mesh into quads, and writes to a .msh file.

    Parameters:
    surfaces (list[int]): List of Gmsh surface tags.
    filename (str): Output .msh file name.
    """

    # Options for generating and saving
    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.option.setNumber("Mesh.SaveParametric", 0)

    # Mesh the surfaces
    for surface in surfaces:
        gmsh.model.geo.mesh.setRecombine(2, surface)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    # Save the mesh to a file
    gmsh.write(filename)


def generate_data(laminar: bool, Re: float, nelq: int = 7, yplus: float = 1.0):
    """
    Generates data to set-up the case. Uses rough yplus-based wall meshing.

    Parameters:
    laminar (bool): is the flow laminar.
    Re (float): what is the flow Reynolds number.
    nelq (int): number of quadrature per element, roughly.
    yplus (float): minimum yplus to be roughly achieved.

    Returns:
    coord_x (list[float]): List of x-coordinates.
    coord_y (list[float]): List of y-coordinates.
    disc_x (list[int]): Discretization values for lines parallel to the x-direction.
    disc_y (list[int]): Discretization values for lines parallel to the y-direction.
    prog_x (list[float]): Progression values for lines parallel to the x-direction.
    prog_y (list[float]): Progression values for lines parallel to the y-direction.
    exclude (list[tuple[int, int]]): List of (i, j) indices to exclude.
    """
    growth_rate = 1.4
    num_y_prog = 4
    num_x = 4

    # With respect to a single arm of the cross
    ymin = -1.0
    ymax = 1.0
    xlen = 2.0

    if laminar:
        ytilde = yplus / np.sqrt(3 * Re)
    else:
        ytilde = yplus * 20 / Re

    y0 = ytilde * nelq

    dyprog = 0
    rate = 1
    for i in range(num_y_prog - 1):
        rate = growth_rate**i
        dyprog += y0 * rate

    rate = growth_rate ** (num_y_prog - 1)
    dy1 = y0 * rate
    dycst = ymax - ymin - 2 * dyprog
    num_y_cst = int(np.ceil(dycst / dy1)) + 1

    # symmetric
    coord_x = [
        -xlen + ymin,
        ymin,
        ymin + dyprog,
        ymin + dyprog + dycst,
        ymax,
        ymax + xlen,
    ]
    coord_y = [
        -xlen + ymin,
        ymin,
        ymin + dyprog,
        ymin + dyprog + dycst,
        ymax,
        ymax + xlen,
    ]

    # exclude corners
    exclude = [
        (0, 0),
        (len(coord_x) - 1, 0),
        (0, len(coord_y) - 1),
        (len(coord_x) - 1, len(coord_y) - 1),
    ]

    # discretization
    disc_x = [num_x, num_y_prog, num_y_cst, num_y_prog, num_x]
    disc_y = [num_x, num_y_prog, num_y_cst, num_y_prog, num_x]

    # element size progression
    prog_x = [1, growth_rate, 1, 1 / growth_rate, 1]
    prog_y = [1, growth_rate, 1, 1 / growth_rate, 1]

    return coord_x, coord_y, disc_x, disc_y, prog_x, prog_y, exclude


# Example usage
if __name__ == "__main__":
    coord_x, coord_y, disc_x, disc_y, prog_x, prog_y, exclude = generate_data(
        LAMINAR, RE
    )

    gmsh.initialize()
    gmsh.model.add("2D Mesh Points and Lines")

    points, point_map = generate_points(coord_x, coord_y, exclude)
    print(f"Generated {len(points)} points: {points}")

    line_map = generate_lines(coord_x, coord_y, point_map, exclude)
    print(f"Generated {len(line_map)} lines.")

    apply_transfinite_curves(
        line_map, coord_x, coord_y, disc_x, disc_y, prog_x, prog_y, exclude
    )
    print(f"Meshed {len(line_map)} lines.")

    surfaces = generate_surfaces(coord_x, coord_y, line_map, exclude)
    print(f"Generated {len(surfaces)} surfaces: {surfaces}")

    add_physical_groups(coord_x, coord_y, line_map, surfaces, exclude)
    print("Added physical groups.")

    mesh_and_save(surfaces, "cross.msh")

    # Finalize Gmsh
    gmsh.finalize()
