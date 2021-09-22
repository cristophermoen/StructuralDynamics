#' # Lecture 5 - Harmonic excitation of dynamic systems

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Harmonic excitation can be a big problem for structural and mechanical systems, especially if $\omega ≈ $ω_n$, which causes dynamic amplification to grow without bounds.  

#' # Let's model harmonic excitation with the forcing function $p(t) = p_osin(\omegat)$.

function equation_of_motion_with_a_forcing_function(utt, ut, u, p, t)
	
    k, m, c, ω, po = p

    utt[1] = -k/m * u[1] - c/m * ut[1] - po/m * sin(ω*t[1])  

end;


#' Apply a harmonic excitation to a W14x38 simply supported steel beam.  Use an SDOF lumped mass and stiffness approximation.

#' Define beam properties.
L = 20.0 #m
E = 1.99948e+11 #N/m^2
I = 385.0 * 25.4^4 /1000^4 #m^4, W14x38 
A = 11.2 * 25.4^2 / 1000^2 #m^2, W14x38
γ = 7700 #kg/m^3, density of steel

#' Calculate the beam lumped mass.
m = γ * A * L  #kg

#' Approximate the beam stiffness.
q = 1.0  #N/m
Δ_midspan = 5 * q * L^4 / (384 * E * I)
P_lumped = q * L  #N   lumped distributed load
k = P_lumped /Δ_midspan  #N/m

#' Define viscous damping parameters.
ωn = sqrt(k/m)   #radians/sec
c_cr = 2 * ωn * m  #kg/sec
ξ = 0.05
c = ξ * c_cr  #kg/sec

#' Define the harmonic forcing function magnitude and frequency.
po = 10000   #N
ω = 0.2 * ωn

#' We don't need an initial velocity or initial displacement to get the beam going since we are using a forcing function.
ut_o = [0.0] 
u_o = [0.0]

#' Package up the important physical parameters.
p = [k, m, c, ω, po];

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 10.0
dt = 0.001

#' Build our model. 
prob = SecondOrderODEProblem(equation_of_motion_with_a_forcing_function, ut_o, u_o, (t_start, t_end), p);

#' And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#' And plot.
using Plots
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "beam disp. u [m]")

#' Calculate the beam natural frequency.
fn = ωn/(2*π)  #cycles/sec.

#' Calculate the forcing function frequency.
f = ω/(2*π)   #cycles/sec.


