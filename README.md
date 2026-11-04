# MATLAB Poisson FDM Solver

This repository contains the public MATLAB script `LAPLACE_FDM.m` for a two-dimensional Poisson problem on a 1 x 1 m square domain.

## Problem definition

The script solves the Poisson equation

```text
 d2u/dx2 + d2u/dy2 = 4
```

with Dirichlet boundary values taken from `u(x,y) = x^2 + y^2` along the square boundary.

## Discretization and linear system

`NX` and `NY`, defined near the beginning of the script, set the number of grid subdivisions along x and y. The defaults are 50 subdivisions in each direction. The script sets `DX = 1/NX` and `DY = 1/NY`, giving `NX + 1` by `NY + 1` nodes.

Second-order finite differences form a five-point stencil for the interior nodes. The resulting coefficient system is represented through five diagonal arrays, `D1` through `D5`, rather than a full matrix.

## Weighted Jacobi iteration

The script prompts for the user-defined relaxation parameter `w` and applies a weighted Jacobi update to the interior values. Iteration continues while the relative algebraic residual is greater than `1e-6`.

## Run in MATLAB

1. Open `LAPLACE_FDM.m` in MATLAB.
2. Run the script.
3. At the grid prompt, enter `1` to plot the computational grid or `0` to skip it.
4. At the relaxation prompt, enter the desired value for `w`.

The script plots the residual history and numerical solution isolines. The number of iterations depends on the selected `w`.

## Scope of validation

This public script does not perform numerical-versus-analytical post-processing or calculate a maximum solution error. The personal portfolio website includes additional post-processing visualizations and fine-grid validation results generated separately. Those website results are not produced by this script.
