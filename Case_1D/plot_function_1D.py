import numpy as np
import matplotlib.pyplot as plt


# Plot a 1D mesh
# Inputs:
    # NodeList: coordinates of the nodes
    # ElemList: connectivity. Each line contains the nodes of an element
def plot_mesh_1D(NodeList,ElemList):

    plt.figure()
    
    # Plot elements
    for i in range(ElemList.shape[0]):
        Nodes = ElemList[i]
    
        for j in range(Nodes.shape[0]-1):
            plt.plot([NodeList[Nodes[j]],NodeList[Nodes[j+1]]],[0.,0.],'-k')

        # Write element number
        plt.text((NodeList[Nodes[0]]+NodeList[Nodes[-1]])/2,0.05,str(i+1),fontweight='bold',ha='center',va='center')
    
    # Plot nodes
    for i in range(NodeList.shape[0]):

        # Difference between inner nodes (within an element) and outer nodes (element boundaries)
        if (i==ElemList.T[0]).any() or (i==ElemList.T[-1]).any():
            alpha = 1
            mec = 'black'
        else:
            alpha = 0.7
            mec = 'none'
            
        plt.plot(NodeList[i,0],0.,'o',markersize=20,color='limegreen',mec=mec,alpha=alpha)
        plt.text(NodeList[i,0],0.,str(i+1),fontweight='bold',ha='center',va='center')
       
    plt.xlim([-0.25,NodeList[:,0].max()+0.25])
    plt.ylim([-0.25,0.25])
    
    return
