using DifferentialEquations 
using LinearAlgebra  

#Import the structural properties of the 3 story mystery building: the stiffness, mass, and damping matrices.
include("MKC.jl")

#' Run some tests here to figure out the natural frequencies, mode shapes, and modal damping ratios.  Modify below as necessary.

function mdof(utt, ut, u, p, t)

  M, K, C = p

  utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:, 1]
  #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end

#' Define the time range for the simulation.
t = 0.0:0.01:5.0

#' Define the problem.
u_o = [0.0; 0.0; 0.0]   
ut_o = [0.0; 0.0; 0.0]

problem = SecondOrderODEProblem(mdof, ut_o, u_o, (0.,5.),(M, K, C))

#' Solve.
solution = solve(problem, DPRKN8(),tstops=t)
