# -*- coding: utf-8 -*-

# Created on Tue Nov  7 19:09:02 2023



import numpy as np
import matplotlib.pyplot as plt
import plot_function_2D as fun

# Entradas.

Lx = 5     # Longitud en eje x [m]
Ly = 5      # Longitud en eje y [m]
nex = 2     # Número de elementos en el eje x
ney = 2     # Número de elementos en el eje y
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

#%% Integración => Kelem

def dEta(xi:float,elem:int)->float:
    if elem == 1:
        deta = (1/4)*(xi-1)
    elif elem==2:
        deta = (-1/4)*(xi+1)
    elif elem==3:
        deta = (1/4)*(xi+1)
    elif elem==4:
        deta = (1/4)*(1-xi)
    
    return deta

def dXi(eta:float,elem:int)->float:
    if elem == 1:
        dxi = (1/4)*(eta-1)
    elif elem==2:
        dxi = (1/4)*(1-eta)
    elif elem==3:
        dxi = (1/4)*(1+eta)
    elif elem==4:
        dxi = (1/4)*(1+eta)*(-1)
    return dxi


def Jacobian(eta:float,xi:float, coordinates:np.ndarray[float])->np.ndarray[float]:
    x = coordinates[:,0]
    y = coordinates[:,1]
    dxdxi = np.matmul(np.array([[dXi(eta,1),dXi(eta,2),dXi(eta,3),dXi(eta,4)]]),x)
    dydxi = np.matmul(np.array([[dXi(eta,1),dXi(eta,2),dXi(eta,3),dXi(eta,4)]]),y)

    dxdeta = np.matmul(np.array([[dEta(xi,1),dEta(xi,2),dEta(xi,3),dEta(xi,4)]]),x)
    dydeta = np.matmul(np.array([[dEta(xi,1),dEta(xi,2),dEta(xi,3),dEta(xi,4)]]),y)

    J = np.array([[dxdxi[0],dydxi[0]],[dxdeta[0],dydeta[0]]])
    
    return J

def MatrixB(eta:float,xi:float, coordinates:np.ndarray[float]):
    J = Jacobian(eta,xi, coordinates)
    J_inv = np.linalg.inv(J)
    dN = np.zeros((4,2))
    for i in range(0,4):
       d = np.transpose(np.matmul(J_inv,np.array([[dXi(eta,i+1)],[dEta(xi,i+1)]])))
       dN[i,:]+=d[0,:]
    
    B = np.array([[dN[0,0],0 , dN[1,0], 0, dN[2,0], 0, dN[3,0], 0], 
                  [0, dN[0,1], 0, dN[1,1], 0, dN[2,1], 0, dN[3,1]],
                  [dN[0,1], dN[0,0], dN[1,1], dN[1,0], dN[2,1], dN[2,0], dN[3,1], dN[3,0]]])
    return B

def MatrixD(E:float,nu:float)->np.ndarray[float]:
    D = E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
    return D

def ComputeElemStiff(coord,D,GPE,t):
    Gp = [[-1/np.sqrt(3),-1/np.sqrt(3)],[1/np.sqrt(3),-1/np.sqrt(3)],[1/np.sqrt(3),1/np.sqrt(3)],[-1/np.sqrt(3),1/np.sqrt(3)]]
    W = 1.0
    Kelem = np.zeros((8,8))
    for i in range(0,GPE):
        xi = Gp[i][0]
        eta = Gp[i][1]
        B = MatrixB(eta,xi,coord)
        Bt = np.transpose(B)
        J = Jacobian(eta,xi,coord)
        detJ = np.linalg.det(J)
        A1 = np.matmul(D,B)
        A2 = np.matmul(Bt,A1)
        Kelem += A2*detJ*t*W
    return Kelem

def testKelem(nodes,elem,D,t):
    elem1 = elem[1,:]
    coord = np.zeros((len(elem1),2))
    for i in range(0,len(elem1)):
        coord[i,:]+=nodes[elem1[i],:]
    kelem = ComputeElemStiff(coord,D,4,t)
    return kelem

#%% Matriz de rigidez global.

# Función para calcular la matriz de rigidez global. 
       
# {INPUTS}: Kglobal, Kelem, elem.
# {OUTPUTS}: Kglobal.

Kglobal = np.zeros((18,18))
Kelem = np.ones((8,8))
elem = elements[0,:]
E = 71e9
nu = 0.33
t = 0.1
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
    for i in range(0,np.shape(elements)[0]):
        print(i)
        elem = np.transpose(elements[i,:])
        Kglobal = Assembly_global_k(Kglobal, Kelem, elem)
    return (Kglobal)
        
SOLVE = iterativeK(Kglobal, Kelem, elements)


test = testKelem(nodes, elements, MatrixD(E,nu), t)