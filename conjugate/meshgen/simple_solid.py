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
laminar = True
nele = 7
Re = 1000
yplus = 1

growth_rate = 1.4
num_y_prog = 5
num_x_ele = 3

ymin = 1.0
ymax = 2.0
xmin = 0.0
xmax = 1.0
zval = 0.0

# -----------------------------------------------------------------------------
# 1) Define geometry
# -----------------------------------------------------------------------------
# Create two points on the x-axis
p1 = gmsh.model.geo.addPoint(xmin, ymin, zval)
p2 = gmsh.model.geo.addPoint(xmin, ymax, zval)

# Create lines for each region
l1 = gmsh.model.geo.addLine(p1, p2)

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
gmsh.model.geo.mesh.setTransfiniteCurve(l1, num_y_prog, "Progression", growth_rate)

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
    [(1, l1)],  # Which entities to extrude (1D line)
    1.0,
    0,
    0,  # (dx, dy, dz) - extrude in +x direction
    numElements=[num_x_ele],  # number of layers in the extrusion
    heights=[xmax],  # total height of the extrusion
    recombine=True,  # ensure quads
)

gmsh.model.geo.synchronize()

# -----------------------------------------------------------------------------
# 4) Add physical groups
# -----------------------------------------------------------------------------
inlet_tags = [l1]
inlet_group = gmsh.model.addPhysicalGroup(1, inlet_tags, 5)
gmsh.model.setPhysicalName(1, inlet_group, "SolidMin")

all_lines = gmsh.model.getEntities(dim=1)  # Get all 1D entities
all_line_tags = [line[1] for line in all_lines]  # Extract tags

# Filter x lines only
x_line_tags = [
    line
    for line in all_line_tags
    if abs(
        gmsh.model.getBoundingBox(1, line)[1] - gmsh.model.getBoundingBox(1, line)[4]
    )
    < 1e-6 * (ymax - ymin)  # y_min ≈ y_max
]

# Dictionary to store x lines and their y-coordinates
line_positions = {}

for line_tag in x_line_tags:
    bbox = gmsh.model.getBoundingBox(1, line_tag)
    print(bbox)
    y_min, y_max = bbox[1], bbox[4]  # Extract y-coordinates
    line_positions[line_tag] = y_min  # Use y_min to sort

print(line_positions)

# Find left and right sides
# Sort by y-coordinate
sorted_lines = sorted(line_positions.items(), key=lambda item: item[1])
left_side = sorted_lines[0][0]  # Line with the smallest y-coordinate
right_side = sorted_lines[-1][0]  # Line with the largest y-coordinate

# Assign to physical groups
left_side_group = gmsh.model.addPhysicalGroup(1, [left_side], 6)
gmsh.model.setPhysicalName(1, left_side_group, "BottomSolid")

right_side_group = gmsh.model.addPhysicalGroup(1, [right_side], 7)
gmsh.model.setPhysicalName(1, right_side_group, "TopSolid")

# The outlet (right) is everything else
outlet_tags = set(all_line_tags) - set(x_line_tags) - set(inlet_tags)
outlet_group = gmsh.model.addPhysicalGroup(1, list(outlet_tags), 8)
gmsh.model.setPhysicalName(1, outlet_group, "SolidMax")

# Add a physical group for the surface
surfaces = [tag[1] for tag in extrusion if tag[0] == 2]  # Extract all surfaces
main_surface_group = gmsh.model.addPhysicalGroup(
    2, surfaces, 102
)  # All extruded surfaces
gmsh.model.setPhysicalName(2, main_surface_group, "SolidSurface")

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
gmsh.write("simple_solid.msh")

# Finalize and exit
gmsh.finalize()
