#' # Lecture 15 - Two story shear building free vibration

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Let's calculate the undamped free vibration response of a two story shear building with Chopra Ch. 10 modal expansion ideas and compare it to a numerical solution.

#' Write the equation of motion.
function mdof(utt, ut, u, p, t)

    M, K = p

    utt[:,1] = - M \ K * u[:,1]
    #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end

#' Describe the building. 

m1 = 10000  #kg   
m2 = 10000  #kg   

k1 = 12000000  #N/m
k2 = 12000000  #N/m 

#' Define mass and stiffness matrices.

M = [m1 0
     0 m2]

K = [k1+k2 -k2
     -k2   k2]

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

#' Show that the scaling factor give {ϕn}'[M]{ϕn} = 1.0

ϕn1 = mode_shapes[:, 1]
ϕn2 = mode_shapes[:, 2]

Mn1 = ϕn1' * M * ϕn1
Mn2 = ϕn2' * M * ϕn2

#' Assume the initial displacement for free vibration is equal to ϕ1 and the initial velocities of both dof are zero.
u_o = [100; 103.6]
ut_o = [0.0; 0.0]

#' Define the time range for the simulation.
t = 0.0:0.01:5.0

#' Calculate the mode 1 contribution first.

#' Calculate the initial conditions.
qn1_o = ϕn1'*M*u_o / (ϕn1'*M*ϕn1)
An1 = qn1_o

qtn1_o = ϕn1'*M*ut_o / (ϕn1'*M*ϕn1)
Bn1 = qtn1_o/ωn[1]

#' Calculate the mode 1 displacements.
un1 = ϕn1' .* (An1 .* cos.(ωn[1].*t) .+ Bn1 .* sin.(ωn[1].*t))

#' Plot.
plot(t, un1)

#' Calculate the mode 2 contribution.

#' Calculate the initial conditions.
qn2_o = ϕn2'*M*u_o / (ϕn2'*M*ϕn2)
An2 = qn2_o

qtn2_o = ϕn2'*M*ut_o / (ϕn2'*M*ϕn2)
Bn2 = qtn2_o/ωn[1]

#' Calculate the mode 1 displacements.
un2 = ϕn2' .* (An2 .* cos.(ωn[2].*t) .+ Bn2 .* sin.(ωn[2].*t))

#' Plot.
plot(t, un2)

#' Define the problem.
#'                                   u_dot0      u_0     trange
problem = SecondOrderODEProblem(mdof, ut_o, u_o, (0.,3.),(M, K, earthquake))

#' Solve.
solution = solve(problem, DPRKN8(),tstops=0:0.01:3.0)

#' Get response.
u1_dot=(x->x[1]).(solution.u)
u2_dot=(x->x[2]).(solution.u)
u1_numerical=(x->x[3]).(solution.u)
u2_numerical=(x->x[4]).(solution.u)
t_numerical=solution.t

#' Plot.
using Plots

plot(t_numerical, u1_numerical)
plot!(t, un1[:, 1] + un2[:, 1])
