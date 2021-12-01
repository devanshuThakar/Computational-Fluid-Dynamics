## Lid driven cavity flow
 A standard benchmark problem of lid-driven cavity flow in computational fluid dynamics is solved. Methods like (Quassi) Rhie-chow interpolation, SIMPLE method and Naive method are coded to solve the problem. A Finite Volume (FV) approch was adopted. 
 
 Rhie-chow interpolation method is used to solve the unsteady Navier-Stokes equation. An explicit euler scheme and second-order central differencing scheme was used for convective terms. A collocated grid for pressure and velocity was used. SIMPLE (Semi-Implicit Method for Pressure Linked Equations) was used to solve the steady state Navier-Stokes equation in a staggered mesh arrangement. 

The velocity field and the boundary conditions are shown in the figure below : 
<p align="center">
    <img src="gifs-images/animation.gif" >
    <!-- <img src="Images/velocity_field.jpg" width="45%"> -->
</p>

<!-- The contour plots of <sub>&phi;</sub> for values of diffusivity constant <sub>&Gamma;</sub>=0.01 and <sub>&Gamma;</sub>=0.001 are shown below : 
<p align="center">
    <img src="Images/QUICK_Gamma_0.001_40X40.png" width="45%">
    <img src="Images/QUICK_Gamma_0.010_40X40.png" width="45%">
</p>

The plots of flux through the west wall with number of CVs by QUICK and CDS scheme is shown below : 
<p align="center">
    <img src="Images/Plot_of_flux_vs_CVs_QUICK.png" width="45%">
    <img src="Images/Plot_of_flux_vs_CVs_CDS.png" width="45%">
</p> -->
