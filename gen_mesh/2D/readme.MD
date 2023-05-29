# 2D Mesh generation

1. `gen_mesh_a.py` : Generate 2D axisymmetric mesh of the nanopippete with triangular elements. It requires the radius $`a_{s}`$ in nm as first argument. The script can be modified for user specific arguments, mesh sizes and dimensions. Outputs gmesh readable `.msh`. Test mesh produced for convergence in PNP in `test_mesh/2D`.

2. `gmsh2xml.py` : Convert `.msh` output to dolfin-xml `.xml` readable by fenics 2019. 