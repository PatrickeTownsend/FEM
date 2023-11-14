import functions as func
#---Variables---#
Lx = 5.0
Ly = 5.0
ne_x = 2
ne_y = 2
t = 0.1
E = 71e9
nu = 0.33
#---Main---#
if __name__ == "__main__":
   elements,nodes = func.discretize(Lx,Ly,ne_x,ne_y)
   #func.plot_mesh(nodes,elements)
   D = func.MatrixD(E,nu)
   kelem = func.testKelem(nodes,elements,D,t)
   print(kelem)

