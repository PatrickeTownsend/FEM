import numpy as np
import matplotlib.pyplot as plt
def discretize(Lx:float,Ly:float,ne_x:int,ne_y:int)->tuple[np.ndarray[int],np.ndarray[float]]:
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
    return np.vstack(elem),np.vstack(nodes)

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

    