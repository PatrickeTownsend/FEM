import numpy as np
import matplotlib.pyplot as plt
def N_nat_fun_1D(xi:np.ndarray[float], NodesPerElem:int) -> np.ndarray[float]:
    """
    Calcula y devuelve las funciones de forma evaluadas para cada valor de xi a lo largo de un elemento en coordenadas naturales

    INPUT
    - `xi` Valores de xi en coordenadas naturales
    - `NodesPerElem` Número de nodos por cada elemento (Determina si es de orden cuadrático o lineal)
    
    OUTPUT
    - `N_matrix` Contiene las funciones de forma evaluadas para cada valor de xi
    """
    # Creo el array vacío #
    N_matrix = np.zeros((NodesPerElem, len(xi)))
    # Determina el orden del elemento (Lineal o cuadratico) #
    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
    # Comienza a evaluar cada funcion de forma y a rellenar al array vacio
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

def discretize(Lx:float,ne:int,NodesPerElem:int)->tuple[np.ndarray[float],np.ndarray[float]]:
    """
    Genera la discretizacion del elemento y devuelve las coordenadas de cada nodo y los nodos de cada elemento

    INPUT
    - `Lx` Longitud del elemento
    - `ne` Número de elementos a discretizar
    - `NodesPerElem` Número de nodos por elemento (Determina si es de orden cuadrático o lineal)

    OUTPUT
    - `nodes` Vector con las coordenadas de cada nodo generado
    - `elements` Array de arrays con los nodos asignados a cada elemento
    """
    # Determina el orden del elemento
    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
    if Linear:
        # determina el numero de nodos total
        nodesNum = ne+1
        # genera el array de nodos
        nodes = np.linspace(0,Lx,nodesNum)
        # crea la matriz de elementos, cada fila es un elemento con un vector que contiene los índices de cada nodo
        elements = np.zeros((ne,2),int)
        for i in range(0,ne):
            # rellena la matriz de elementos
            elements[i,:]+=[i,i+1]
    elif cuadratic:
        nodesNum = 2*ne+1
        nodes = np.linspace(0,Lx,nodesNum)
        elements = np.zeros((ne,3),int)
        for i in range(0,ne):
            elements[i,:]=[2*i,2*i+1,2*i+2]
    return nodes,elements


def plot_mesh_1D(NodeList:np.ndarray[float], ElemList:np.ndarray[float]):
    """
    Grafica el mallado con la discretizacion creada

    INPUT
    - `nodeList` Vector con las coordenadas de cada nodo generado
    - `ElemList` Array de arrays con los nodos asignados a cada elemento
    
    """
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

def interpolation(N_matrix:np.ndarray[float],nodes:np.ndarray[float],elements:np.ndarray[float],NodesPerElem:int)->tuple[np.ndarray[float],np.ndarray[float]]:
    """
    Interpola los desplazamiento en cada nodo de la malla generada

    INPUT
    - `N_matrix` Contiene las funciones de forma evaluadas para cada valor de xi
    - `nodes` Vector con las coordenadas de cada nodo generado
    - `elements` Array de arrays con los nodos asignados a cada elemento
    - `NodesPerElem` Número de nodos por elemento 

    OUTPUT
    - `x_global` Vector con las coordenadas globales a lo largo de todo el dominio
    - `u_global` Vector con los desplazamientos globales interpolados a lo largo de todo el dominio
    """
    # Determina el orden de cada elemento
    if NodesPerElem == 2:
        Linear = True
        cuadratic = False
    else:
        cuadratic = True
        Linear = False
    # Inicializa los arrays
    NumElement = np.shape(elements)[0] # Obtiene el numero de elementos creados
    N_matrix = np.transpose(N_matrix) #Adapta la matriz a la dimension deseada [Nx2] - N: longitud del vector xi
    x_global = np.zeros((NumElement,np.shape(N_matrix)[0])) # Matriz de [ne,N]
    u_global = np.zeros((NumElement,np.shape(N_matrix)[0])) # Matriz de [ne,N]
    for i in range(0,NumElement):
        start = elements[i][0] # primer nodo del elemento
        end = elements[i][1] # segundo nodo del elementos
        x1 = nodes[start] # coordenada del primer nodo
        x2 = nodes[end] # coordenada del segundo nodo
        u1 = Desp_function(x1) # desplazamiento del primer nodo
        u2 = Desp_function(x2) # desplazamiento del segundo nodo
        if Linear:
           x_local = np.array([[x1],[x2]]) # genera el vector de [2x1]
           u_local = np.array([[u1],[u2]]) # genera un vector de [2x1]
        elif cuadratic:
           x3 = nodes[elements[i][2]] # coordenada del tercer nodo (si es cuadratico)
           u3 = Desp_function(x3) # desplazamiento del tercer nodo (si es cuadratico)
           x_local = np.array([[x1],[x2],[x3]]) # genera el vector de [3x1]
           u_local = np.array([[u1],[u2],[u3]]) # genera el vector de [3x1]
        x = np.dot(N_matrix,x_local) # Interpolacion de coordenadas [N][x]
        u = np.dot(N_matrix,u_local) # Interpolacion de desplazamientos [N][u]
        for j in range(0,np.shape(N_matrix)[0]):
            x_global[i,j]+=x[j]
            u_global[i,j]+=u[j]
        
    x_global = x_global.flatten() # genera el vector interpolado
    u_global = u_global.flatten() # genera el vector interpolado
    return x_global, u_global

def Desp_function(x:float)->float:
    """
    Calcula el desplazamiento exacto para un punto x del dominio

    INPUT
    - `x` coordenada del nodo

    OUTPUT
    - `u` desplazamiento del punto x
    """
    u = -4*(x-0.3)**6 + (0.6-x)**5 - (5*10**(-6))*(1/(x+0.1))**4 + (x-0.2)**3
    return u

def u_exact(x:np.ndarray[float])->np.ndarray[float]:
    """
    Genera un vector con el despalzamiento exacto a lo largo de todo el dominio

    INPUT
    - `x` Nube de puntos a lo largo de todo el dominio

    OUTPUT
    - `u` Desplazamientos a lo largo de todo el dominio
    """
    u = np.zeros((len(x)))
    for i in range(0,len(x)):
        u[i] += Desp_function(x[i]) # evalua la funcion en cada punto de x
    return u

def errorDist(u_global:np.ndarray[float],Lx:float)->tuple[np.ndarray[float],np.ndarray[float]]:
    """
    Genera la distribucion del error de interpolacion sobre todo el dominio

    INPUT
    - `u_global` Vector con los desplazamientos globales interpolados a lo largo de todo el dominio
    - `Lx` Longitud todo el dominio

    OUTPUT
    - `error` error de interpolacion en cada punto del dominio
    - `x` malla de puntos a lo largo de todo el dominio
    """
    x = np.linspace(0,Lx,len(u_global))
    u_e = u_exact(x) # calcula los desplazamientos exactos para todo el dominio
    error = np.zeros(len(u_global)) 
    for i in range(0,len(u_global)):
        error[i] = (abs(u_global[i]-u_e[i])/abs(u_e[i])) # error en cada punto de malla
    return error,x

def errorRelTol(u_global:np.ndarray[float],Lx:float)->float:
    """
    Calcula el error relativo de interpolacion para todo el dominio

    INPUT
    - `u_global` Vector con los desplazamientos globales interpolados a lo largo de todo el dominio
    - `Lx` Longitud total del dominio

    OUTPUT
    - `errorAbs` Error de interpolacion de todo el dominio
    """
    x = np.linspace(0,Lx,len(u_global))
    u_e = u_exact(x)
    errorAbs = (np.linalg.norm(u_global - u_e)/np.linalg.norm(u_e))
    return errorAbs

def errorRefinement(xi:np.ndarray[float],Lx:float,NodesPerElem:int,typeOfElem:str):
    """
    Itera sobre varios valores de ne (numero de elementos) y grafica el error relativo en funcion de los GDL

    INPUT
    - `xi` Valores de xi en coordenadas naturales
    - `Lx` Longitud total de todo el dominio
    - `NodesPerElem` Número de nodos por cada elemento (Determina si es de orden cuadrático o lineal)
    - `typeOfElem` String con el tipo de elemento (lineal o cuadrático)
    """
    ne = 1
    errorAbs = 1.0
    epsilon=[]
    DOF =[]
    GDL = 0
    while errorAbs >= 0.02: # Itera para cada valor de ne y fija la tolerancia deseada
       N_matrix = N_nat_fun_1D(xi,NodesPerElem)
       nodes,elements = discretize(Lx,ne,NodesPerElem)
       x_global,u_global = interpolation(N_matrix,nodes,elements,NodesPerElem)
       errorAbs = errorRelTol(u_global,Lx)
       epsilon.append(errorAbs)
       # Computa el número de GDL total para todo el dominio (actualiza para valor de ne)
       if ne ==1:
           if NodesPerElem==2:
               GDL += 2
           else:
               GDL += 3
       else:
           if NodesPerElem==2:
               GDL+=1
           else:
               GDL+=2
       DOF.append(GDL)
       ne+=1
    # Imprime el numero de GDL y elementos para un error por debajo del 2%
    print("---------------------------------------------------------------------------")
    print("Se necesitan ",ne-1," elementos ",typeOfElem,"con ",GDL,"GDL para obtener un error de ",errorAbs)
    print("---------------------------------------------------------------------------")  
    # Grafica el error vs GDL
    plt.figure()
    plt.plot(DOF,epsilon)
    plt.axhline(y=0.02,color='r',label="error del 2%")
    title = "Error relativo vs GDL elementos "+ typeOfElem
    plt.title(title)
    plt.xlabel("Números de GDLs")
    plt.ylabel("Error relativo")
    plt.grid()
    plt.legend()
    route = "plots/error_vs_gdl_"+typeOfElem+".png"
    plt.savefig(route)
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
    route="plots/interpolation_"+typeOfElem+".png"
    plt.savefig(route)
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
    route = "plots/form_func_"+typeOfElem+".png"
    plt.savefig(route)
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