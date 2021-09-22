
#' # Lecture 7 - Simulating a blast loading 

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Blast loadings are impulses with typically some small back pressure (suction).

#' Let's consider a typical blast pressure record on a facade wall of a building.  I received this record from Dr. Matthew Whelan at the University of North Carolina - Charlotte.  

#' The pressure record units are psi, and the time interval dt = 0.0001 seconds.

using DelimitedFiles
blast_pressure = vec(readdlm("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_7/blast_pressure.txt"))  
t = collect(range(0,length=5001,stop=0.5))  

using Plots
plot(t, blast_pressure)
xlabel!("time [sec.]")
ylabel!("blast pressure [psi]")

#' The blast pressure record spikes above 40 psi and then you can see a small backpressure (-0.1 psi) as the blast shock wave recedes.

#' Model facade wall behavior as an SDOF dynamic system.  

#' Write the equation of motion.

function equation_of_motion(utt, ut, u, p, t)
	
    k, m, c, blast_loading_function = p

    utt[1] = -k/m * u[1] - c/m * ut[1] - blast_loading_function(t[1])

end;

#' Assume a wall tributary width.
tributary_width = 1.0  #m

#' Assume a wall height.
wall_height = 3.0 #m

#' Calculate the effective wall area.
A = tributary_width * wall_height  #m^2

#' Convert the blast pressure from psi to N/m^2.
blass_pressure = blast_pressure * 6894.76

#' Calculate the blast force from the pressure.
blast_force = blast_pressure .* A  #N

#' Define the blast loading as a continuous function.
using Dierckx
blast_loading_function = Spline1D(t, blast_force)

#' Assume a mass for the effective wall area.
m = 2500 #kg

#' Assume a linear stiffness for now.  The stiffness is most likely nonlinear as the vertical framing member transitions from flexure to catenary action.   
k = 4448/0.050  #N/m

#' Define viscous damping parameters.
ωn = sqrt(k/m)   #radians/sec
c_cr = 2 * ωn * m  #kg/sec
ξ = 0.05
c = ξ * c_cr  #kg/sec

#' We don't need an initial velocity or initial displacement to get the wall going since we are using a forcing function.
ut_o = [0.0] 
u_o = [0.0]

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 10.0
dt = 0.001

#' Package up the important physical parameters.
p = [k, m, c, blast_loading_function];

#' Build our model. 
prob = SecondOrderODEProblem(equation_of_motion, ut_o, u_o, (t_start, t_end), p);

#' And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#' And plot.
using Plots
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [m]")