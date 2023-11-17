import functions as func
import numpy as np
import matplotlib.pyplot as plt
import plot_functions_2D_deformed as plotf
#---Variables---#
Lx = 5.0
Ly = 5.0
ne_x = 50
ne_y = 50
t = 0.1
E = 71e9
nu = 0.33
side_cond = "left" # (right,left,upper,lower)
side_force = "upper"
DOF = 1 # (1 if vertical) / (0 if horizontal)
delta = 0
F_dist = 1e10
F_pun = 0
#---Main---#
if __name__ == "__main__":
   elements,nodes,mesh = func.discretize(Lx,Ly,ne_x,ne_y)
   D = func.MatrixD(E,nu)

   #---Generating arrays-----#
   Kglobal = func.StiffnesMatrix(nodes,elements,D,t)
   Fglobal = np.zeros((np.shape(Kglobal)[0],1))

   #----Penalty----#
   Boundary_cond = func.LocateNodes("lower",mesh,1)
   Boundary_cond2 = func.LocateNodes("lower",mesh,0)
   Boundary_cond3 = func.LocateNodes("left",mesh,0)
   Boundary_force = func.LocateNodes("upper",mesh,1)
   #print(Kglobal[0:4,0:4])

   #----Boundary conditions----#
   Kb = np.max(Kglobal)*10e6
   Kglobal,Fglobal = func.ApplyDirichletPen(Boundary_cond,delta,Kb,Kglobal,Fglobal)
   Kglobal,Fglobal = func.ApplyDirichletPen(Boundary_cond2,delta,Kb,Kglobal,Fglobal)
   Kglobal,Fglobal = func.ApplyDirichletPen(Boundary_cond3,delta,Kb,Kglobal,Fglobal)
   #print(Kglobal[0:4,0:4])

   #----Penalty-----#
   Fglobal= func.ApplyNewman(Boundary_force,F_pun,Fglobal)
   Fglobal = func.FuerzaDist(Boundary_force,F_dist,Fglobal,t,Lx)
   #print(Fglobal)

   #----Solving----#
   U = np.linalg.solve(Kglobal,Fglobal)
   #print(U)
   nodes_new = func.Update_nodes(nodes,U)
   func.plot_mesh(nodes,elements)
   plotf.plot_deformed_mesh(nodes,elements,nodes_new,Plot_nodes=True,Plot_elem=True)
   U_new = U.reshape(U.shape[0])
   plotf.PlotDisplacements(nodes,elements,U_new,Plot_elem=False)

  


   

   

