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
        elements = np.zeros((ne,2),int)
        for i in range(0,ne):
            elements[i,:]+=[i,i+1]
    elif cuadratic:
        nodesNum = 2*ne+1
        nodes = np.linspace(0,Lx,nodesNum)
        elements = np.zeros((ne,3),int)
        for i in range(0,ne):
            elements[i,:]=[2*i,2*i+1,2*i+2]
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
    plt.title("Discretización del dominio")

    return

def interpolation(N_matrix,nodes,elements,NodesPerElem):
    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
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
        if Linear:
           x_local = np.array([[x1],[x2]])
           u_local = np.array([[u1],[u2]])
        elif cuadratic:
           x3 = nodes[elements[i][2]]
           u3 = Desp_function(x3)
           x_local = np.array([[x1],[x2],[x3]])
           u_local = np.array([[u1],[u2],[u3]])
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

def errorDist(u_global,Lx):
    x = np.linspace(0,Lx,len(u_global))
    u_e = u_exact(x)
    error = np.zeros(len(u_global))
    for i in range(0,len(u_global)):
        error[i] = (abs(u_global[i]-u_e[i])/abs(u_e[i]))
    return error,x

def errorRelTol(u_global,Lx):
    x = np.linspace(0,Lx,len(u_global))
    u_e = u_exact(x)
    errorAbs = (np.linalg.norm(u_global - u_e)/np.linalg.norm(u_e))
    return errorAbs

def errorRefinement(xi,Lx,NodesPerElem,typeOfElem):
    ne = 2
    errorAbs = 1.0
    epsilon=[]
    GDL = []
    while errorAbs > 0.02:
       N_matrix = N_nat_fun_1D(xi,NodesPerElem)
       nodes,elements = discretize(Lx,ne,NodesPerElem)
       x_global,u_global = interpolation(N_matrix,nodes,elements,NodesPerElem)
       errorAbs = errorRelTol(u_global,Lx)
       epsilon.append(errorAbs)
       GDL.append(len(nodes)+1)
       ne += 1
    print("---------------------------------------------------------------------------")
    print("Se necesitan ",ne," elementos ",typeOfElem,"con ",len(nodes)+1,"GDL para obtener un error de ",errorAbs)
    print("---------------------------------------------------------------------------")  
    plt.figure()
    plt.plot(GDL,epsilon)
    plt.axhline(y=0.02,color='r',label="error del 2%")
    title = "Error relativo vs GDL elementos "+ typeOfElem
    plt.title(title)
    plt.xlabel("Número de GDL")
    plt.ylabel("Error relativo")
    plt.grid()
    plt.legend()
    return

def Plot_Interpolated(x,x_global,u_global,nodes,typeOfElem):
    plt.figure()
    plt.plot(x,u_exact(x),label="Exacta",linestyle="--")
    plt.plot(x_global,u_global,label = "Interpolada")
    plt.plot(nodes,np.zeros(len(nodes)),marker="o")
    plt.scatter(nodes,u_exact(nodes),color="r",marker="o",zorder=3,label="Desplazamiento nodal")
    title = "Desplazamientos elementos " + typeOfElem
    plt.title(title)
    plt.xlabel("Posicion (m)")
    plt.ylabel("Desplazamiento")
    plt.legend()
    plt.grid()
    return

def Plot_formFunc(xi,N_matrix,NodesPerElem,typeOfElem):
    plt.figure()
    plt.plot(xi,N_matrix[0,:],label="Funcion N1")
    plt.plot(xi,N_matrix[1,:],label="Funcion N2")
    if NodesPerElem ==3:
       plt.plot(xi,N_matrix[2,:],label="Funcion N3")
    title = "Funciones de Forma elemento "+ typeOfElem
    plt.title(title)
    plt.legend()
    plt.xlabel("xi")
    plt.grid()
    return

def PlotErrorDist(u_global,Lx):
    error,x_e = errorDist(u_global,Lx)
    plt.figure()
    plt.plot(x_e,error)
    plt.axhline(y=0.02,color='r')
    plt.title("Distribucion de error")
    plt.xlabel("Posicion (m)")
    plt.ylabel("error relativo")
    plt.grid()