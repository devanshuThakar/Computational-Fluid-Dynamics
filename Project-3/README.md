## Scalar transport by Finite Volume
A 2D scalar transport equation for the transport of quantiy <sub>&phi;</sub> without source term and at steady state conditions was solved using Finite Volume method (FV). The boundary conditions and the velocity field are given. The three popular aprroches of discretizing the convectine term was compared. They are linear interpolation scheme (CDS), upwind difference scheme (UDS) and the QUICK scheme. 

The velocity field and the boundary conditions are shown in the figure below : 
<p align="center">
    <img src="Images/Grids.jpg" width="45%">
    <!-- <img src="Images/velocity_field.jpg" width="45%"> -->
</p>

The contour plots of <sub>&phi;</sub> for values of diffusivity constant <sub>&Gamma;</sub>=0.01 and <sub>&Gamma;</sub>=0.001 are shown below : 
<p align="center">
    <img src="Images/QUICK_Gamma_0.001_40X40.png" width="45%">
    <img src="Images/QUICK_Gamma_0.010_40X40.png" width="45%">
</p>

The plots of flux through the west wall with number of CVs by QUICK and CDS scheme is shown below : 
<p align="center">
    <img src="Images/Plot_of_flux_vs_CVs_QUICK.png" width="45%">
    <img src="Images/Plot_of_flux_vs_CVs_CDS.png" width="45%">
</p>
