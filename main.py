import numpy as np
import matplotlib.pyplot as plt
import functions as func
#----Input Variables---#
NodesPerElem = 3
Lx = 1.0
ne = 5
xi = np.linspace(-1,1,100)
x = np.linspace(0,Lx,100)
N_matrix = func.N_nat_fun_1D(xi,NodesPerElem)
nodes,elements = func.discretize(Lx,ne,NodesPerElem)
print(elements)
x_global,u_global = func.interpolation(N_matrix,nodes,elements)
plt.figure()
plt.plot(x,func.u_exact(x),label="Exacta")
plt.plot(x_global,u_global,label = "Interpolada")
plt.plot(nodes,np.zeros(len(nodes)),marker="o")
plt.legend()
#----Plots----#
func.plot_mesh_1D(nodes,elements)
plt.figure()
plt.plot(xi,N_matrix[0,:],label="Funcion N1")
plt.plot(xi,N_matrix[1,:],label="Funcion N2")
if NodesPerElem ==3:
    plt.plot(xi,N_matrix[2,:],label="Funcion N3")
plt.title("Funciones de Forma")
plt.legend()
plt.show()
