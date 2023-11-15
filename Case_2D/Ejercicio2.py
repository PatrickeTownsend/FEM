# -*- coding: utf-8 -*-

# Created on Tue Nov  7 19:09:02 2023



import numpy as np
import matplotlib.pyplot as plt
import plot_function_2D as fun

# Entradas.

Lx = 5     # Longitud en eje x [m]
Ly = 5      # Longitud en eje y [m]
nex = 2     # Número de elementos en el eje x
ney = 1     # Número de elementos en el eje y
t = 0.1     # Espesor [m]
F = 10e6    # Carga distribuida [N/m] (+: hacia la superficie)


#%% Discretización del dominio.

# Función de discretización: discretiza un dominio dimensional en nodos
# y elementos. Tipo de elementos aceptado: cuadriláteros bilineales
# de tensión plana.

# {INPUTS}: Lx, Ly, nex, ney.
# {OUTPUTS}: nodes, elements.
    
def discretize(Lx, Ly, nex, ney):
    deltax = Lx / nex                    # Distancia entre nodos del eje x
    deltay = Ly / ney                    # Distancia entre nodos del eje y
    nodes = []                           # Inicializa nodes
    elements = []                        # Inicializa elements
      
    # Se van rellenando las coordenadas de los nodos
    for j in range(ney + 1):
      for i in range(nex + 1):
          x = i*deltax
          y = j*deltay
          nodes.append([x, y])
    nodes = np.array(nodes)  
    
    # Se forma la matriz de conectividad
    for j in range(ney):
      for i in range(nex):
          n1 = i+(nex+1)*j   	# Columna 1
          n2 = n1+1             # Columna 2
          n3 = n1+nex+1         # Columna 3
          n4 = n3+1             # Columna 4
          elements.append([n1, n2, n4, n3])   # Genera la matriz
    elements = np.array(elements)
    
    return nodes, elements

# Representación de nodos y elementos.

nodes, elements = discretize(Lx, Ly, nex, ney)
fun.plot_mesh(nodes,elements)

#%% Matriz de rigidez global.

# Función para calcular la matriz de rigidez global. 
       
# {INPUTS}: Kglobal, Kelem, elem.
# {OUTPUTS}: Kglobal.

Kglobal = np.zeros((12,12))
Kelem = np.ones((8,8))
elem = elements[0,:]
# for i in range(0,8):
#     for j in range(0,8):
#         Kelem[i,j] = j+1+(1+i)*10
# elem = elements[1,:]


def Assembly_global_k(Kglobal, Kelem, elem):
    Nnodes = len(elem)
    for i in range(Nnodes):
        for j in range(Nnodes):
            Kglobal[2*elem[i]:2*(elem[i]+1), 2*elem[j]:2*(elem[j]+1)] += Kelem[i*2:i*2+2, j*2:j*2+2]
    return(Kglobal)


def iterativeK(Kglobal,Kelem,elements):
    for i in range(0,len(elements)):
        elem = elements[i,:]
        Kglobal += Assembly_global_k(Kglobal, Kelem, elem)
    return (Kglobal)
        
SOLVE = Assembly_global_k(Kglobal, Kelem, elem)
    
    
    