# 2DNSSolver_SpectralDealiased

## Description

A spectral Legendre collocation method using pseudo time
stepping that solves the steady Navier-Stokes equations
for the lid driven cavity flow problem (sub-critical
Reynolds numbers).

For marginally resolved flows error from aliasing becomes a problem. 
The code use de-aliasing of the convective terms through a projection approach,
which is helpful for making an otherwise diverging solution, converge.

## Example results

Solution on 48x48 grid: Contours of stream function (left) and contours of vorticity
with indication of magnitude; negative to positive colored blue to red (right).
![Stream function and vorticity contours](https://github.com/spietz/2DNSSolver_SpectralDealiased/blob/master/fig/ContoursStreamVort.png)

Comparison of extreme velocities at centerlines
![Comparison with Ghia et al.](https://github.com/spietz/2DNSSolver_SpectralDealiased/blob/master/fig/ComparisonGhia.png)

## References

D.A.Kopriva, Implementing Spectral Methods for partial differential
equations, Springer, 2009.  

U. Ghia, K.N. Ghia, C.T. Shin, High-Re solutions for incompressible 
flow using the Navier-Stokes equations and a multigrid method, 
J. Comput. Phys. 1982.

W. Zhang, C.H. Xi, An explicit Chebyshev pseudospectral multigrid method for 
incompressibe Navier-Stokes equations, Computers & Fluids, 2010.
