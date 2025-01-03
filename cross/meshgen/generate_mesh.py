import gmsh

# Initialize Gmsh
gmsh.initialize()
# Optional: display messages in the terminal
gmsh.option.setNumber("General.Terminal", 1)

# Add a new model
gmsh.model.add("non_uniform_line_extrude")

# -----------------------------------------------------------------------------
# 1) Define geometry
# -----------------------------------------------------------------------------
# Create two points on the x-axis
p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
p2 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0)

# Create a line between p1 and p2
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

num_segments = 20
growth_rate = 1.2  # For example, larger near p2
gmsh.model.geo.mesh.setTransfiniteCurve(l1, num_segments, "Progression", growth_rate)

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
    [(1, l1)],  # Which entity to extrude (1D line)
    0,
    1.0,
    0,  # (dx, dy, dz) - extrude in +y direction
    numElements=[10],  # number of layers in the extrusion
    heights=[1.0],  # total height of the extrusion
    recombine=True,  # ensure quads
)

gmsh.model.geo.synchronize()

# -----------------------------------------------------------------------------
# 4) Enforce quads-only (recombine all surfaces)
# -----------------------------------------------------------------------------
gmsh.option.setNumber("Mesh.RecombineAll", 1)
# Optional: The transfinite setup also helps ensure structured/quadrilateral meshes

# -----------------------------------------------------------------------------
# Generate and save
# -----------------------------------------------------------------------------
gmsh.model.mesh.generate(2)
gmsh.write("non_uniform_extruded_quads.msh")

# Finalize and exit
gmsh.finalize()
