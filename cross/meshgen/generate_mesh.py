import gmsh

# Initialize Gmsh
gmsh.initialize()
# Optional: display messages in the terminal
gmsh.option.setNumber("General.Terminal", 1)

# Add a new model
gmsh.model.add("non_uniform_line_extrude")

# -----------------------------------------------------------------------------
# 0) Data
# -----------------------------------------------------------------------------
xmin = 0.0
xmax = 1.0
ymin = 0.0
ymax = 1.0
zval = 0.0

dxprog = 0.4

growth_rate = 1.2
num_x = 20
num_y = 10

# -----------------------------------------------------------------------------
# 1) Define geometry
# -----------------------------------------------------------------------------
# Create two points on the x-axis
p1 = gmsh.model.geo.addPoint(xmin, ymin, zval)
p2 = gmsh.model.geo.addPoint(xmax, ymin, zval)

# Define additional points to split the line into three segments
p_mid1 = gmsh.model.geo.addPoint(
    dxprog, ymin, zval
)  # Transition point to constant region
p_mid2 = gmsh.model.geo.addPoint(
    xmax - dxprog, ymin, zval
)  # Transition point to decreasing region

# Create lines for each region
l1 = gmsh.model.geo.addLine(p1, p_mid1)  # Increasing mesh size
l2 = gmsh.model.geo.addLine(p_mid1, p_mid2)  # Constant mesh size
l3 = gmsh.model.geo.addLine(p_mid2, p2)  # Decreasing mesh size

# -----------------------------------------------------------------------------
# 2) Impose a non-uniform 1D mesh distribution
# -----------------------------------------------------------------------------
#
# setTransfiniteCurve(lineTag, numberOfNodes, progression)
#
# - numberOfNodes: total number of divisions along the line
# - progression:   how the element size varies along the line
#                  > 1   -> smaller elements near the start point
#                  < 1   -> smaller elements near the end point
#                  = 1.0 -> uniform distribution

# Set transfinite curve for increasing mesh size
gmsh.model.geo.mesh.setTransfiniteCurve(l1, num_x, "Progression", growth_rate)

# Set transfinite curve for constant mesh size
gmsh.model.geo.mesh.setTransfiniteCurve(l2, num_x, "Progression", 1.0)

# Set transfinite curve for decreasing mesh size (mirror of l1)
gmsh.model.geo.mesh.setTransfiniteCurve(l3, num_x, "Progression", 1.0 / growth_rate)

# Make the geometry entities consistent
gmsh.model.geo.synchronize()

# -----------------------------------------------------------------------------
# 3) Extrude the line into 2D (mimicking a wall-bounded style mesh)
# -----------------------------------------------------------------------------
# extrude(dimTagged, dx, dy, dz, [numElements], [heights], recombine)
#
# - dimTagged = (1, l1) refers to a 1D entity with tag = l1
# - Here we extrude in the y-direction by 1.0, with 10 divisions
# - recombine=True ensures quadrilateral elements in the newly created surface

extrusion = gmsh.model.geo.extrude(
    [(1, l1), (1, l2), (1, l3)],  # Which entities to extrude (1D line)
    0,
    1.0,
    0,  # (dx, dy, dz) - extrude in +y direction
    numElements=[num_y],  # number of layers in the extrusion
    heights=[ymax],  # total height of the extrusion
    recombine=True,  # ensure quads
)

gmsh.model.geo.synchronize()

# -----------------------------------------------------------------------------
# 4) Add physical groups
# -----------------------------------------------------------------------------
inlet_tags = [l1, l2, l3]
inlet_group = gmsh.model.addPhysicalGroup(1, inlet_tags)
gmsh.model.setPhysicalName(1, inlet_group, "Inlet")

all_lines = gmsh.model.getEntities(dim=1)  # Get all 1D entities
all_line_tags = [line[1] for line in all_lines]  # Extract tags

# Filter vertical lines only
vertical_line_tags = [
    line
    for line in all_line_tags
    if abs(
        gmsh.model.getBoundingBox(1, line)[0] - gmsh.model.getBoundingBox(1, line)[3]
    )
    < 1e-6 * (xmax - xmin)  # x_min ≈ x_max
]

# Dictionary to store vertical lines and their x-coordinates
line_positions = {}

for line_tag in vertical_line_tags:
    bbox = gmsh.model.getBoundingBox(1, line_tag)
    x_min, x_max = bbox[0], bbox[3]  # Extract x-coordinates
    line_positions[line_tag] = x_min  # Use x_min to sort

# Find left and right sides
# Sort by x-coordinate
sorted_lines = sorted(line_positions.items(), key=lambda item: item[1])
left_side = sorted_lines[0][0]  # Line with the smallest x-coordinate
right_side = sorted_lines[-1][0]  # Line with the largest x-coordinate

# Assign to physical groups
left_side_group = gmsh.model.addPhysicalGroup(1, [left_side])
gmsh.model.setPhysicalName(1, left_side_group, "LeftSide")

right_side_group = gmsh.model.addPhysicalGroup(1, [right_side])
gmsh.model.setPhysicalName(1, right_side_group, "RightSide")

# The outlet (top) is everything else
outlet_tags = set(all_line_tags) - set(vertical_line_tags) - set(inlet_tags)
outlet_group = gmsh.model.addPhysicalGroup(1, list(outlet_tags))
gmsh.model.setPhysicalName(1, outlet_group, "Outlet")

# Add a physical group for the surface
surfaces = [tag[1] for tag in extrusion if tag[0] == 2]  # Extract all surfaces
main_surface_group = gmsh.model.addPhysicalGroup(2, surfaces)  # All extruded surfaces
gmsh.model.setPhysicalName(2, main_surface_group, "MainSurface")

# -----------------------------------------------------------------------------
# 5) Set options for saving
# -----------------------------------------------------------------------------
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.ElementOrder", 2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.option.setNumber("Mesh.Binary", 0)
gmsh.option.setNumber("Mesh.SaveAll", 0)
gmsh.option.setNumber("Mesh.SaveParametric", 0)

# -----------------------------------------------------------------------------
# Generate and save
# -----------------------------------------------------------------------------
gmsh.model.mesh.generate(2)
gmsh.write("non_uniform_extruded_quads.msh")

# Finalize and exit
gmsh.finalize()
