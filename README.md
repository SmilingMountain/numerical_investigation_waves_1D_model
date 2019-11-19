# numerical_investigation_waves_1D_model

## A set of a Matlab scripts to numerically solve the equations for the behaviour of waves in a one-dimensional shallow water model, under periodic and open boundary conditions

This project investigates the behaviour of waves in a one-dimensional shallow water model by using
a series of MATLAB scripts to numerically solve the equations involved. It starts with a linearised version of the
shallow water differential equations and analyses the evolution of the system both mechanically and
energetically under periodic and open boundary conditions. 

This is realized by implementing an
Euler forward step and a "leapfrog" scheme in the MATLAB code. In general, differentials can be 
transformed to centred differences, forward differences and backward
differences, as presented in the full report. The Euler forward step uses a forward difference in time 
and a centred difference in space, while the "leapfrog" forward step uses centred differences 
in both space and time. 

The complete, non-linear system, is then considered and again analysed 
under the two sets of boundary conditions, including the more complex case of orography.

The shallow water equations, which are derived from the Navier-Stokes equations under the assumption
that the horizontal length scale is much larger than the vertical one, are a set of hyperbolic
partial differential equations that describe the evolution of a hydrostatic homogeneous incompressible 
fluid in response to gravitational and rotational accelerations. The one-dimensional equations 
are used extensively in computer models because they are significantly easier to solve than
the full shallow water equations.

The full project report contains all the necessary background information and details, methods used, results and conclusions.
The MATLAB script files contain code for the various cases considered, as indicated by their names, and include an appropriate level of commenting and documentation.


