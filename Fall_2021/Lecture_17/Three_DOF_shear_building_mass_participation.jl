#' # Lecture 17 - Three story shear building spatial mass distribution and participation

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Calculate the spatial distribution of mass in each mode.

#' Describe the building. 

m1 = 10000  #kg   
m2 = 10000  #kg   
m3 = 10000  #kg

k1 = 12000000  #N/m
k2 = 12000000  #N/m 
k3 = 12000000  #N/m 

#' Define mass and stiffness matrices.

M = [m1 0 0
     0 m2 0
     0 0 m3]

K = [k1+k2 -k2 0
     -k2   k2+k3  -k3
     0     -k3   k3]

#' Calculate the mode shapes.
using LinearAlgebra

#' Solve for eigenvalues of K*ϕn=wn^2*M*ϕn.
ωn_squared=eigvals(K, M)

#' Calculate the modal natural frequencies.
ωn = sqrt.(ωn_squared)

#' Calculate the mode shape natural frequencies in Hz.
fn = ωn ./(2 * π)

#' Calculate the mode shape natural periods in seconds.
Tn = 1 ./ fn

#' Solve for eigenvectors K*ϕn=wn^2*M*ϕn.
mode_shapes=eigvecs(K,M)

ϕ1 = mode_shapes[:, 1]
ϕ2 = mode_shapes[:, 2]
ϕ3 = mode_shapes[:, 3]


#' Calculate the spatial mass distribution for each of the modes.

#' Define the mass influence vector.  Since all degrees of freedom are lateral, then for an earthquake, then the mass at each story will be involved.  

ι = [1.0, 1.0, 1.0]  

#' Calculate the modal mass spatial (per story) distribution vector, and the scale factor for this vector (per mode).  
M * ϕ1
ϕ1' * M * ι

M * ϕ2
ϕ2' * M * ι

M * ϕ3
ϕ3' * M * ι


sn1 = ϕ1' * M * ι * M * ϕ1

sn2 = ϕ2' * M * ι * M * ϕ2

sn3 = ϕ3' * M * ι * M * ϕ3

s = sn1 + sn2 + sn3

#' Calculate the modal mass participation ratio.

sum(sn1)

ρ1 = sum(sn1)/(m1+m2+m3)

ρ2 = sum(sn2)/(m1+m2+m3)

ρ3 = sum(sn3)/(m1+m2+m3)

ρ = ρ1 + ρ2 + ρ3