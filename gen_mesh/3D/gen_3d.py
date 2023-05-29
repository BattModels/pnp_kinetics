import numpy as np
from mpi4py import MPI
import warnings
warnings.filterwarnings("ignore")
import gmsh, os, sys

gmsh.initialize()
gmsh.model.add("Model1")
#gmsh.model.setCurrent("Cylinder")
#channel = gmsh.model.occ.addBox(0, 0, 0, L, B, H)
a_s = 1e-7

Z = 2e-4
r = 1e-6
R = 3.5e-5
print("Z={}, r={}, R={}".format(Z,r,R))
cylinder = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, r, r)
frustum = gmsh.model.occ.addCone(0, 0, r, 0, 0, Z-r, r, R)
fluid = gmsh.model.occ.fuse([(3, cylinder)], [(3, frustum)])
gmsh.model.occ.synchronize()

## Mark volume
volumes = gmsh.model.getEntities(dim=3)
assert(volumes == fluid[0])
fluid_marker = 11
gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], fluid_marker)
gmsh.model.setPhysicalName(volumes[0][0], fluid_marker, "Fluid volume")

## Dilate mesh to factor
factor = 0.1 # Transform mesh dimensions to factor
gmsh.model.occ.dilate(gmsh.model.getEntities(), 0 ,0 , 0, factor, factor, factor)
Z = 2e-4 * factor
r = 1e-6 * factor
R = 3.5e-5 * factor

## Assign a mesh size to all the points:
#gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lcar1)

# Mark surfaces
surfaces = gmsh.model.occ.getEntities(dim=2)
#bot_id, top_id, cyl_wall_id, cone_wall_id = 1, 3, 5, 7
for surface in surfaces:
    com = gmsh.model.occ.getCenterOfMass(surface[0], surface[1])
    #if np.allclose(com, [0, 0, 0], rtol=0.1):
    if np.allclose(com[2], 0, rtol=0.1):
        #gmsh.model.addPhysicalGroup(surface[0], [surface[1]], bot_id)
        bot = surface[1]
        #gmsh.model.setPhysicalName(surface[0], bot_id, "Bottom Surface")
    elif np.allclose(com[2], r/2, rtol=0.1):
        #gmsh.model.addPhysicalGroup(surface[0], [surface[1]], cyl_wall_id)
        cyl_wall = surface[1]
        #gmsh.model.setPhysicalName(surface[0], cyl_wall_id, "Cylinder Wall")
    elif np.allclose(com[2], Z, rtol=0.1):
        #gmsh.model.addPhysicalGroup(surface[0], [surface[1]], top_id)
        top = surface[1]
        #gmsh.model.setPhysicalName(surface[0], top_id, "Top Surface")
    else:
        #gmsh.model.addPhysicalGroup(surface[0], [surface[1]], cone_wall_id)
        cone_wall = surface[1]
        #gmsh.model.setPhysicalName(surface[0], cone_wall_id, "Cone Wall")

gmsh.model.mesh.field.add("Distance", 1)
gmsh.model.mesh.field.setNumbers(1, "FacesList", [top]) # Mesh bw bottom and top wall
gmsh.model.mesh.field.setNumbers(1, "FacesList", [bot])
gmsh.model.mesh.field.setNumber(1, "Sampling", 1000)
res = 1*r #1*r orig
gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(2, "IField", 1)
gmsh.model.mesh.field.setNumber(2, "LcMin", 0.24*res) #0.3*res orig
gmsh.model.mesh.field.setNumber(2, "LcMax", 10*res) #10*res orig
gmsh.model.mesh.field.setNumber(2, "DistMin", 1*r)
gmsh.model.mesh.field.setNumber(2, "DistMax", Z)

## Additional mesh size specifications
# dist2 = gmsh.model.mesh.field.add("Distance")
# gmsh.model.mesh.field.setNumbers(dist2, "FacesList", [cone_wall])
# res = 10000*r
# thres2 = gmsh.model.mesh.field.add("Threshold")
# gmsh.model.mesh.field.setNumber(thres2, "IField", dist2)
# gmsh.model.mesh.field.setNumber(thres2, "LcMin", res)
# gmsh.model.mesh.field.setNumber(thres2, "LcMax", 4*res)
# gmsh.model.mesh.field.setNumber(thres2, "DistMin", 0.2*res)
# gmsh.model.mesh.field.setNumber(thres2, "DistMax", 0.5*res)

minimum = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", [2])#, thres2])
gmsh.model.mesh.field.setAsBackgroundMesh(minimum)

def meshSizeCallback(dim, tag, x, y, z, lc):
    return min(lc, 0.02 * x + 0.01)

#p1 = gmsh.model.occ.addPoint(0, 0, 0)
#p2 = gmsh.model.occ.addPoint(0, 1e-7, 0)
#p3 = gmsh.model.occ.addPoint(0.7e-7, 0.69e-7, 0)
gmsh.model.occ.synchronize()

gmsh.model.mesh.generate(3)
gmsh.write("mesh3D.msh")

os.system('meshio convert mesh3D.msh mesh3D.xml')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
