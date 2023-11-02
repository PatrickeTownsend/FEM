import numpy as np
import matplotlib.pyplot as plt
import functions as func
node_number = 3
xi = np.linspace(-1,1,1000)

N_matrix = func.N_nat_fun_1D(xi,node_number)


plt.plot(xi,N_matrix[0,:],label="Funcion N1")
plt.plot(xi,N_matrix[1,:],label="Funcion N2")
if node_number ==3:
    plt.plot(xi,N_matrix[2,:],label="Funcion N3")
plt.title("Funciones de Forma")
plt.legend()
plt.show()
