import gmsh
import sys
import numpy as np

A = int(sys.argv[1]) # Radius of nano-pipette in nm

gmsh.initialize()

gmsh.model.add("mesh_{}".format(str(int(A))))

a = A*1e-9
l = 2e-5 # Length of pipette
rs = (a-1e-7) + 3.5e-6 # Radius at the top

def inter_point(point1, point2, frac=0.5):
	return [frac*(point2[0]-point1[0]), frac*(point2[1]-point1[1])]


int_p = inter_point([a, a], [rs, l], 1/3)

lc0 = 8e-8
lc1 = (1e-2)*1e-8
lc2 = 2e-8

print(lc0, lc1, lc2)

p1 = gmsh.model.geo.addPoint(0, 0, 0, lc1)
p2 = gmsh.model.geo.addPoint(a, 0, 0, lc1)
p3 = gmsh.model.geo.addPoint(a, a, 0, lc1)
p4 = gmsh.model.geo.addPoint(int_p[0], int_p[1], 0, lc2)
p5 = gmsh.model.geo.addPoint(rs, l, 0, lc0)
p6 = gmsh.model.geo.addPoint(0, l, 0, lc0)
p7 = gmsh.model.geo.addPoint(0, a, 0, lc1)

l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p5)
l5 = gmsh.model.geo.addLine(p5, p6)
l6 = gmsh.model.geo.addLine(p6, p7)
l7 = gmsh.model.geo.addLine(p7, p1)

#cl1 = gmsh.model.geo.addCurveLoop([p1, p2, p3, p4, p5])
#s1 = gmsh.model.geo.addPlaneSurface([cl1])

gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6, l7], 1)
gmsh.model.geo.addPlaneSurface([1], 1)


gmsh.model.geo.synchronize()

#gmsh.model.addPhysicalGroup(1, [p1, p2, p4], pg1)
#gmsh.model.addPhysicalGroup(2, [s1], name = "My surface")

gmsh.model.mesh.generate(2)

# save it to disk
gmsh.write("mesh_"+str(int(A))+"_test.msh")


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
