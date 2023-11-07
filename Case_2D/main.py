import functions as func
#---Variables---#
Lx = 5.0
Ly = 5.0
ne_x = 5
ne_y = 3
#---Main---#
if __name__ == "__main__":
   elements,nodes = func.discretize(Lx,Ly,ne_x,ne_y)
   func.plot_mesh(nodes,elements)
