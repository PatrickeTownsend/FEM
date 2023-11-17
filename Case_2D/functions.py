import numpy as np
import matplotlib.pyplot as plt
def discretize(Lx:float,Ly:float,ne_x:int,ne_y:int)->tuple[np.ndarray[int],np.ndarray[float],np.ndarray[float]]:
    """
    Genera la malla del elementos discretizado

    INPUT
    - `Lx` Longitud en x del dominio
    - `Ly` Longitud en y del dominio
    - `ne_x` número de particiones en x
    - `ne_y` número de particiones en y

    OUTPUT
    - `nodes` Array con las coordenadas de cada nodo
    - `elements` Array con los nodos de cada elemento
    """
    nodes = []
    elem = []
    mesh = np.zeros((ne_y+1,ne_x+1),int)
    dx = Lx/ne_x
    dy = Ly/ne_y
    node=0
    for j in range(0,ne_y+1):
        for i in range(0,ne_x+1):
            mesh[j,i] += node
            nodes.append([i*dx,j*dy])
            node += 1
    for j in range(0,ne_y):
        for i in range(0,ne_x):
            elem_mesh = mesh[j:j+2,i:i+2]
            elem.append([elem_mesh[0,0],elem_mesh[0,1],elem_mesh[1,1],elem_mesh[1,0]])
    return np.vstack(elem),np.vstack(nodes),mesh

# Plot a 2D mesh
# Inputs:
    # nodes: coordinates of the nodes. Each row is a node, with its x and y coordinates in the first and second column, respectively.
    # elements: connectivity matrix. Each row contains the node indices of an element.
def plot_mesh(nodes,elements):

    plt.figure()
    
    # Plot elements
    for i in range(elements.shape[0]):
        Nodes = elements[i]
        Node1 = nodes[int(Nodes[0])]
        Node2 = nodes[int(Nodes[1])]
        Node3 = nodes[int(Nodes[2])]
        Node4 = nodes[int(Nodes[3])]
    
        plt.plot([Node1[0],Node2[0]],[Node1[1],Node2[1]],'-k')
        plt.plot([Node2[0],Node3[0]],[Node2[1],Node3[1]],'-k')
        plt.plot([Node3[0],Node4[0]],[Node3[1],Node4[1]],'-k')
        plt.plot([Node4[0],Node1[0]],[Node4[1],Node1[1]],'-k')
    
        plt.text((Node1[0]+Node2[0])/2,(Node1[1]+Node4[1])/2,str(i+1),fontweight='bold',ha='center',va='center')
    
    # Plot nodes
    for i in range(nodes.shape[0]):
        plt.plot(nodes[i,0],nodes[i,1],'o',markersize=20,color='limegreen')   
        plt.text(nodes[i,0],nodes[i,1],str(i+1),fontweight='bold',ha='center',va='center')
            
    plt.xlim([-0.5,np.max(nodes[:,0])+0.5])
    plt.ylim([-0.5,np.max(nodes[:,1])+0.5])

    plt.show()

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

def MatrixB(eta:float,xi:float, coordinates:np.ndarray[float])->np.ndarray[float]:
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

def ComputeElemStiff(coord:np.ndarray[float],D:np.ndarray[float],GPE:int,t:float)->np.ndarray[float]:
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

def ElemCoords(elem:np.ndarray[float],nodes:np.ndarray[float])->np.ndarray[float]:
    coord = np.zeros((len(elem),2))
    for i in range(0,len(elem)):
        coord[i,:]+=nodes[elem[i],:]
    return coord

def Asembly_global(Kglobal:np.ndarray[float],Kelem:np.ndarray[float],elem:np.ndarray[float])->np.ndarray[float]:
    for i in range(len(elem)):
        for j in range(len(elem)):
            Kglobal[2*elem[i]:2*(elem[i]+1), 2*elem[j]:2*(elem[j]+1)] += Kelem[2*i:i*2+2, 2*j:2*j+2]
    return Kglobal

def StiffnesMatrix(nodes:np.ndarray[float],elements:np.ndarray[float],D:np.ndarray[float],t:float)->np.ndarray[float]:
    node_num = np.shape(nodes)[0]
    Kglobal = np.zeros((node_num*2,node_num*2))
    for i in range(np.shape(elements)[0]):
        elem = elements[i,:]
        coord = ElemCoords(elem,nodes)
        Kelem = ComputeElemStiff(coord,D,4,t)
        Kglobal = Asembly_global(Kglobal,Kelem,elem)
    return Kglobal

def LocateNodes(side:str,mesh:np.ndarray[int],DOF:int)->list:
    size = np.shape(mesh)
    Boundary = []
    if side=="right":
        line = mesh[:,size[1]-1]
    elif side=="lower":
        line = mesh[0,:]
    elif side=="left":
        line = mesh[:,0]
    elif side=="upper":
        line = mesh[size[0]-1,:]
    for i in range(len(line)):
        Boundary.append([line[i],DOF])
    return Boundary


def ApplyDirichletPen(Boundary,delta,Kb,Kglobal,Fglobal):
    K_new = np.copy(Kglobal)
    F_new = np.copy(Fglobal)
    for i in range(len(Boundary)):
        index = Boundary[i][0]
        DOF = Boundary[i][1]
        K_new[2*index+DOF,2*index+DOF] += Kb
        F_new[2*index+DOF] = Kb*delta
    return K_new,F_new

def ApplyNewman(Boundary,F,Fglobal):
    F_newman = np.copy(Fglobal)
    for i in range(len(Boundary)):
        index = Boundary[i][0]
        DOF = Boundary[i][1]
        F_newman[2*index+DOF] += F
    return F_newman

def FuerzaDist(Boundary,F,F_newman):
    w = F/(len(Boundary)-1)
    for i in range(len(Boundary)):
        index = Boundary[i][0]
        DOF = Boundary[i][1]
        if i==0 or i==len(Boundary)-1:
            F_newman[2*index+DOF] += w/2
        else:
            F_newman[2*index+DOF] += w
    return F_newman

def Update_nodes(nodes,U):
    desp = U.reshape(np.shape(nodes))
    nodes_desp = nodes + desp
    return nodes_desp







    







    




  

