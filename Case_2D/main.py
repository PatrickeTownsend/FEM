import functions as func
import numpy as np
import matplotlib.pyplot as plt

#---Variables---#
Lx = 10.0
Ly = 1.0
ne_x = 5
ne_y = 1
t = 0.1
E = 71e9
nu = 0.33
side_cond = "left" # (right,left,upper,lower)
side_force = "upper"
DOF = 1 # (1 if vertical) / (0 if horizontal)
delta = 0
F_dist = -10e6
F_pun = 0
#---Main---#
if __name__ == "__main__":
   elements,nodes,mesh = func.discretize(Lx,Ly,ne_x,ne_y)
   func.plot_mesh(nodes,elements)
   D = func.MatrixD(E,nu)

   #---Generating arrays-----#
   Kglobal = func.StiffnesMatrix(nodes,elements,D,t)
   Fglobal = np.zeros((np.shape(Kglobal)[0],1))

   #----Penalty----#
   Boundary_cond = func.LocateNodes(side_cond,mesh,1)
   Boundary_cond2 = func.LocateNodes(side_cond,mesh,0)
   Boundary_force = func.LocateNodes(side_force,mesh,1)

   #----Boundary conditions----#
   Kb = np.max(Kglobal)*10e6
   Kglobal,Fglobal = func.ApplyDirichletPen(Boundary_cond,delta,Kb,Kglobal,Fglobal)
   K_new,F_dir = func.ApplyDirichletPen(Boundary_cond2,delta,Kb,Kglobal,Fglobal)
   Fglobal= func.ApplyNewman(Boundary_force,F_pun,Fglobal)
   F = func.FuerzaDist(Boundary_force,F_dist,Fglobal)

   #----Solving----#
   U = np.linalg.solve(K_new,F)
   nodes_new = func.Update_nodes(nodes,U)
  


   

   

