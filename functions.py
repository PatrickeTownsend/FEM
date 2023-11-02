import numpy as np
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
