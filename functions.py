import numpy as np
import matplotlib.pyplot as plt
def N_nat_fun_1D(xi, NodesPerElem):

    N_matrix = np.zeros((NodesPerElem, len(xi)))
    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
    for i in range(0,NodesPerElem):
        for j in range(0,len(xi)):
            if Linear:
                if i==0:
                    N_matrix[i,j]+= 0.5*(1-xi[j])
                else:
                    N_matrix[i,j]+=0.5*(1+xi[j])
            if cuadratic:
                if i==0:
                    N_matrix[i,j]+=(-xi[j]/2)*(1-xi[j])
                if i==1:
                    N_matrix[i,j]+=(1-xi[j])*(1+xi[j])
                if i==2:
                    N_matrix[i,j]+=(xi[j]/2)*(1+xi[j])

    return N_matrix

def discretize(Lx,ne,NodesPerElem):

    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
    if Linear:
        nodesNum = ne+1
        nodes = np.linspace(0,Lx,nodesNum)
        elements = np.zeros((nodesNum-1,2),int)
        for i in range(0,nodesNum-1):
            elements[i,:]+=[i,i+1]
    elif cuadratic:
        nodesNum = 2*ne+1
        nodes = np.linspace(0,Lx,nodesNum)
        elements = np.zeros((nodesNum-1,2),int)
        for i in range(0,nodesNum-1):
            elements[i,:]+=[i,i+1]
    return nodes,elements


def plot_mesh_1D(NodeList, ElemList):
    plt.figure()

    # Plot elements
    for i in range(ElemList.shape[0]):
        Nodes = ElemList[i]

        for j in range(Nodes.shape[0] - 1):
            plt.plot([NodeList[Nodes[j]], NodeList[Nodes[j + 1]]], [0., 0.], '-k')

        # Write element number
        plt.text((NodeList[Nodes[0]] + NodeList[Nodes[-1]]) / 2, 0.05, str(i + 1), fontweight='bold', ha='center',
                 va='center')

    # Plot nodes
    for i in range(NodeList.shape[0]):

        # Difference between inner nodes (within an element) and outer nodes (element boundaries)
        if (i == ElemList.T[0]).any() or (i == ElemList.T[-1]).any():
            alpha = 1
            mec = 'black'
        else:
            alpha = 0.7
            mec = 'none'

        plt.plot(NodeList[i], 0., 'o', markersize=20, color='limegreen', mec=mec, alpha=alpha)
        plt.text(NodeList[i], 0., str(i + 1), fontweight='bold', ha='center', va='center')

    plt.xlim([-0.25, NodeList[:].max() + 0.25])
    plt.ylim([-0.25, 0.25])

    return

def interpolation(N_matrix,nodes,elements):
    NumElement = np.shape(elements)[0]
    N_matrix = np.transpose(N_matrix)
    x_global = np.zeros((NumElement,np.shape(N_matrix)[0]))
    u_global = np.zeros((NumElement,np.shape(N_matrix)[0]))
    for i in range(0,NumElement):
        start = elements[i][0]
        end = elements[i][1]
        x1 = nodes[start]
        x2 = nodes[end]
        u1 = Desp_function(x1)
        u2 = Desp_function(x2)
        x_local = np.array([[x1],[x2]])
        u_local = np.array([[u1],[u2]])
        x = np.dot(N_matrix,x_local)
        u = np.dot(N_matrix,u_local)
        for j in range(0,np.shape(N_matrix)[0]):
            x_global[i,j]+=x[j]
            u_global[i,j]+=u[j]
        
    x_global = x_global.flatten()
    u_global = u_global.flatten()
    return x_global, u_global

def Desp_function(x):
    u = -4*(x-0.3)**6 + (0.6-x)**5 - (5*10**(-6))*(1/(x+0.1))**4 + (x-0.2)**3
    return u

def u_exact(x):
    u = np.zeros((len(x)))
    for i in range(0,len(x)):
        u[i] += Desp_function(x[i])
    return u
