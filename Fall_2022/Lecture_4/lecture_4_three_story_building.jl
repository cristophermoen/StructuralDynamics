#' # Lecture 4 - Three story shear building free vibration
#' EN 560.630 Fall 2022
#' Department of Civil and Systems Engineering
#' Johns Hopkins University

using DifferentialEquations

#' Let's calculate the undamped free vibration response of a three story shear building with Raleigh damping.

#' Write the equation of motion.
function mdof(utt, ut, u, p, t)

    M, K, C = p

    utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:, 1]
    #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end

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

#' Calculate Raleigh damping C = ao*M + a1*K by assuming the modal damping ratio ζn for modes 1 and 3.

ζ1 = 0.05
ζ3 = 0.05

#' Find ao and a1.
a = (1/2 * [1/ωn[1] ωn[1]
            1/ωn[3] ωn[3]])  \ [ζ1; ζ3] 

#' Define the Raleigh damping matrix.
C = a[1] * M + a[2] * K

#' Let's calculate the modal damping ratio in mode 2, just to make sure things make sense.

ϕ2 = mode_shapes[:, 2]

C2 = ϕ2' * C * ϕ2

M2 = ϕ2' * M * ϕ2

Ccr2 = 2 * ωn[2] * M2

ζ2 = C2 / Ccr2

#' Assume the initial displacement for free vibration is equal to ϕ2 and the initial velocities of all dof are zero.
u_o = ϕ2    
ut_o = [0.0; 0.0; 0.0]

#' Define the time range for the simulation.
t = 0.0:0.01:5.0

#' Define the problem.
#'                                   u_dot0      u_0     trange
problem = SecondOrderODEProblem(mdof, ut_o, u_o, (0.,5.),(M, K, C))

#' Solve.
solution = solve(problem, DPRKN8(),tstops=t)

#' Plot solution.
using Plots
plot(solution)

#' Get response.
u1_dot=(x->x[1]).(solution.u)
u2_dot=(x->x[2]).(solution.u)
u3_dot=(x->x[3]).(solution.u)
u1_numerical=(x->x[4]).(solution.u)
u2_numerical=(x->x[5]).(solution.u)
u3_numerical=(x->x[6]).(solution.u)
t_numerical=solution.t

#' Plot.
using Plots

plot(t_numerical, u1_numerical)
plot!(t_numerical, u2_numerical)
plot!(t_numerical, u3_numerical)