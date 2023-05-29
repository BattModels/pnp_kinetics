import numpy as np
import meshio
import sys

# For 2D mesh only

f_name = sys.argv[1]

geometry = meshio.read(f_name + ".msh")
#meshio.write("mesh.xdmf", meshio.Mesh(points=geometry.points, cells={"triangle": geometry.cell_data_dict['triangle']}))

#cd=MeshFunction('size_t',mesh,loc+'mesh_physical_region.xml')
#fd=MeshFunction('size_t',mesh,loc+'mesh_facet_region.xml')

def createMesh(mesh, cell_type, prune_z=True):
    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:geometrical", cell_type)
    points = mesh.points[:,:2] if prune_z else mesh.points
    out_mesh = meshio.Mesh(points=points, cells={cell_type: cells}, cell_data={"name_to_read":[cell_data]})
    return out_mesh

#line_mesh = createMesh(geometry, "line", prune_z=True)
#meshio.write("facet_test.xdmf", line_mesh)

triangle_mesh = createMesh(geometry, "triangle", prune_z=True)
meshio.write(f_name+".xml", triangle_mesh)
