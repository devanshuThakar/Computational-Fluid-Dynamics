## Lid driven cavity flow
 A standard benchmark problem of lid-driven cavity flow in computational fluid dynamics is solved. Methods like (Quassi) Rhie-chow interpolation, SIMPLE method and Naive method are coded to solve the problem. A Finite Volume (FV) approch was adopted. 
 
 Rhie-chow interpolation method is used to solve the unsteady Navier-Stokes equation. An explicit euler scheme and second-order central differencing scheme was used for convective terms. A collocated grid for pressure and velocity was used. SIMPLE (Semi-Implicit Method for Pressure Linked Equations) was used to solve the steady state Navier-Stokes equation in a staggered mesh arrangement. 

The velocity field and the boundary conditions are shown in the figure below : 
<p align="center">
    <img src="gifs-images/animation.gif" >
</p>