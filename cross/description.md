## Cross
In this case, the goal is to simulate a domain with two inlets and outlets in a cross shape. The mesh is multi-block and built with Gmsh via the Python API. Mesh generation is done in `meshgen`, the other folders correspond to simulation folders.

The main case files are marked as `cross` (e.g. `cross.py` for mesh generation) but two additional, simpler cases are included, `simple_rectangle` and `nonuniform_rectangle`. These have been done first as tests of the Gmsh/Nek5000 workflow.

No attempt is done here to describe the mesh generation procedure in detail. For the two `rectangle` cases, the strategy is to mesh lines and then to generate the mesh via extrusion. For the `cross` cases, a more advanced approach is used. In fact, the resulting set of methods in `cross.py` can be used to generate 2D multi-block meshes with high degree of generality. The symmetric cross example used here is only a single configuration made possible by this approach.

Attention is directed especially to the options for mesh generation and saving. They are necessary for making the resulting mesh Nek5000-compliant.The actual workflow is to generate the `.msh` files with the Python script and then running the Nek5000 utility `gmsh2nek` with the queries answered as:
```
2
<name of case without .msh>
0
0
<name of case without .msh>
```
The resulting `.re2` is to be copied to the corresponding simulation folder, `genmap` can then be used to generate the `.ma2` file.

The simulation case files (`SIZE`, `.par`, `.usr`) do not feature many unique aspects different from typical Nek5000 simulations aside from the `usrdat` approach to setting boundary conditions at the boundaries with their IDs corresponding to those assigned in the mesh generation. Furthermore, the `cross_turbulent` case implements RANS turbulence modelling as per the Nek5000 tutorial with the added modification of inlet boundary conditions for the added RANS scalar quantities.