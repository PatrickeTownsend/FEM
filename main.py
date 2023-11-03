import numpy as np
import matplotlib.pyplot as plt
import functions as func

#----Input Variables---#
N = 200 # Tamaño de la malla 
NodesPerElem = 2 # Number of nodes per element (2: lineal y 3: cuadratico)
ne = 5 # Numero de elementos a discretizar

Lx = 1.0 # Longitud total
xi = np.linspace(-1,1,N) # vector de xi
x = np.linspace(0,Lx,N) # vector x

#---Main Execution---#
if __name__ == "__main__":
   if NodesPerElem==3:
      typeOfElem = "Cuadraticos"
   else:
      typeOfElem = "Lineales"
   N_matrix = func.N_nat_fun_1D(xi,NodesPerElem)
   nodes,elements = func.discretize(Lx,ne,NodesPerElem)
   x_global,u_global = func.interpolation(N_matrix,nodes,elements,NodesPerElem)
   func.errorRefinement(xi,Lx,NodesPerElem,typeOfElem)

   #---Plots---#
   func.Plot_Interpolated(x,x_global,u_global,nodes,typeOfElem)
   func.plot_mesh_1D(nodes,elements)
   func.Plot_formFunc(xi,N_matrix,NodesPerElem,typeOfElem)
   func.PlotErrorDist(u_global,Lx)
   plt.show()

