import matplotlib.pyplot as plt
import numpy as np

# Plot a 2D mesh: before and after displacements
# Inputs:
    # nodes: coordinates of the nodes. Each row is a node, with its x and y coordinates in the first and second column, respectively.
    # elements: connectivity matrix. Each row contains the node indices of an element.
    # nodes_new: updated coordinates of the nodes (after displacement). Each row is a node, with its x and y coordinates in the first and second column, respectively.

def plot_deformed_mesh(nodes,elements,nodes_new,Plot_nodes=True,Plot_elem=True):

    fig_mesh, ax_mesh = plt.subplots(layout='constrained')
    
    # Original mesh
    if Plot_elem:
        Plot_elements(ax_mesh,nodes,elements,linewidth = 0.5, alpha = 0.5)
        
    if Plot_nodes:
        for i in range(nodes.shape[0]):
            ax_mesh.plot(nodes[i,0],nodes[i,1],'o',markersize=5,color='limegreen',alpha=0.5)
    
    # Deformed mesh
    if Plot_elem:
        Plot_elements(ax_mesh,nodes_new,elements,linewidth = 1., alpha = 1.)

    if Plot_nodes:
        for i in range(nodes_new.shape[0]):
            ax_mesh.plot(nodes_new[i,0],nodes_new[i,1],'o',markersize=5,color='red',alpha=1.)
    plt.show()
    



# Plot displacement fields
# Inputs:
    # nodes: coordinates of the nodes. Each row is a node, with its x and y coordinates in the first and second column, respectively.
    # elements: connectivity matrix. Each row contains the node indices of an element.
    # u_sol: displacements. Solution of the FEM system

def PlotDisplacements(nodes,elements,u_sol,Plot_elem=False):

    xv = np.zeros([2,2])
    yv = np.zeros([2,2])
    zv_ux = np.zeros([2,2])
    zv_uy = np.zeros([2,2])

    # Get u_x aund u_y from all nodes
    u_x_sol = u_sol.reshape(nodes.shape[0],2)[:,0]
    u_y_sol = u_sol.reshape(nodes.shape[0],2)[:,1]

    levels_u_x = np.linspace(u_x_sol.min(), u_x_sol.max(), 20)
    levels_u_y = np.linspace(u_y_sol.min(), u_y_sol.max(), 20)
    
    fig_ux, ax_ux = plt.subplots(layout='constrained')
    fig_uy, ax_uy = plt.subplots(layout='constrained')

    for elem in range(elements.shape[0]):
        
        nodes_elem = nodes[elements[elem]]
        
        xv[0] = nodes_elem[0:2,0]
        xv[1] = np.flip(nodes_elem[2:4,0])

        yv[0] = nodes_elem[0:2,1]
        yv[1] = np.flip(nodes_elem[2:4,1])        

        ux = u_sol[2*elements[elem]]
        uy = u_sol[2*elements[elem]+1]

        zv_ux[0] = ux[0:2]
        zv_ux[1] = np.flip(ux[2:4])

        zv_uy[0] = uy[0:2]
        zv_uy[1] = np.flip(uy[2:4])
       
        CS_ux = ax_ux.contourf(xv,yv,zv_ux,levels_u_x,cmap = 'jet')
        CS_uy = ax_uy.contourf(xv,yv,zv_uy,levels_u_y,cmap = 'jet')

    ax_ux.set_title('Displacement u_x')
    ax_uy.set_title('Displacement u_y')
    fig_ux.colorbar(CS_ux)
    fig_uy.colorbar(CS_uy)
    
    
    # Deformed mesh
    if Plot_elem:
        Plot_elements(ax_ux,nodes,elements,linewidth = 1., alpha = 1.)
        Plot_elements(ax_uy,nodes,elements,linewidth = 1., alpha = 1.)
    plt.show()


    return
    
    
# Plot element lines in the given axis ax.
# To be used with plot_deformed_mesh and PlotDisplacements
def Plot_elements(ax,nodes,elements,linewidth = 1., alpha = 1.):

    for i in range(elements.shape[0]):
        Nodes = elements[i]
        Node1 = nodes[int(Nodes[0])]
        Node2 = nodes[int(Nodes[1])]
        Node3 = nodes[int(Nodes[2])]
        Node4 = nodes[int(Nodes[3])]
        
        ax.plot([Node1[0],Node2[0]],[Node1[1],Node2[1]],'-k',linewidth=linewidth,alpha=alpha)
        ax.plot([Node2[0],Node3[0]],[Node2[1],Node3[1]],'-k',linewidth=linewidth,alpha=alpha)
        ax.plot([Node3[0],Node4[0]],[Node3[1],Node4[1]],'-k',linewidth=linewidth,alpha=alpha)
        ax.plot([Node4[0],Node1[0]],[Node4[1],Node1[1]],'-k',linewidth=linewidth,alpha=alpha)
       
    return
