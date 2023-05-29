# 3D Mesh generation

1. `gen_3d.py` : Generate 3D nanopipette mesh with tetrahedral elements. Due tolerance issues, gmesh is first instructed to generate mesh with 10 times bigger dimensions. The mesh is transformed (dilated) to the required dimension later. Mesh tolerance, dimensions, size distribution can be tweaked as per user convenience. Test mesh produced for convergence in PNP in `test_mesh/3D`.  

2. For gmesh conversion to dolfin-xml use e.g. : `meshio convert mesh3D.msh mesh3D.xml`

## Notes
* Do not use gms2xml.py from `2D/`
* Do not have elements other than tetrahedrons in the 3D `.msh` file, `.xml` does not support it