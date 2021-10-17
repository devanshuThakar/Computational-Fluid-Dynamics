## 2D Convection Diffusion problem
A 2D convection diffusion equation with soruce term was solved in this problem. This project is divided into two parts : part-A and part-B.

Part-A solves the steady-state 2D convection diffusion equation with a given source term and boundary conditions. 
![alt-txt](Images/2D-ss.jpg)


The steady-state equation was solved with wide varity of numeircal methods like Gauss Elimination, Jacobi method, Gauss Seidel method, row-wise and column-wise sweep and ADI method. A few representatiove plots from the project are shown below : 
<table>
  <tr>
    <td><img src="Images/CPU_time_vs_grid_size.png" width=270 height=480></td>
    <td><img src="Images/row_column_ADI_residual_Iteration_PartA_Q4.png" width=270 height=480></td>
  </tr>
 </table>
<!-- ![alt-text-1](Images/CPU_time_vs_grid_size.png) | ![alt-text-2](row_column_ADI_residual_Iteration_PartA_Q4.png) -->

Part-B solves the unsteat-state 2D convection diffusion equation with a given source term and boundary conditions. 
![alt-txt](Images/2D-uns.jpg)


The unsteady equation was solved using various numerical methods like the explicit method, implicit method and Crank Nicolson method. 

![alt-text-1](Contour_implicit_t=10_PartB_Q3.png) 

For deatails see the report.pdf file. 