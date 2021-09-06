#' # Lecture 3 - Introduction to damping in dynamic systems

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Let's explore damping in structures. As [Chopra](https://www.pearson.com/us/higher-education/program/Chopra-Dynamics-of-Structures-5th-Edition/PGM1101746.html) says, there can be many different energy dissipating mechanisms that act simulataneously.  There are few popular ones: viscous damping, friction damping, and material damping.  Let's explore them.


#' # Viscous damping

function equation_of_motion_with_viscous_damping(utt, ut, u, p, t)
	
    k, m, c = p

    utt[1] = -k/m * u[1] - c/m * ut[1]  

end;


#' Consider the physical model of the [CFS-NEES building](https://www.youtube.com/watch?v=RXBi66k-Xo0).  The total mass of the building is about [78000 lbf](https://www.ce.jhu.edu/cfsnees/publications/RR01_CFS-NEES%20RR01%20Building%20Structural%20Design%20Narrative.pdf).
m = 78000 * 0.453592; #kg

#' Assume that when you push laterally on the building with 1000 lbf, that it displaces 3.0 in.
k = (1000 * 4448)/(3.0 * 25.4/1000) ; #N/m

#' Model the building as a single degree of freedom in free vibration.  We pull the building laterally over by 15 mm.
ut_o = [0.0]
u_o = [15.0/1000] #mm

#' Assume a viscous damping coefficient by scaling the mass.  The coefficient $c$ has units of kg/sec.
c = 100.0 * m

#' Package up the important physical parameters.
p = [k, m, c];

#' Build our model. 
prob = SecondOrderODEProblem(equation_of_motion_with_viscous_damping, ut_o, u_o, (0.0, 1.0), p);

#' And solve.  
solution = solve(prob, DPRKN6(), tstops=0:0.001:1.0);

#' And plot.
using Plots
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "lateral disp. u [m]")

#' Plot the damping force in the equation of motion as a function of lateral building displacement.
ut = (x->x[1]).(solution.u) #velocity
fd = ut .* c
plot(u, fd, legend = false, xlabel="lateral disp. u [m]", ylabel = "fd [N]")


#' # Friction damping

function equation_of_motion_with_friction_damping(utt, ut, u, p, t)
	
    k, m, fd_friction = p

    if ut[1] >= 0.0

        utt[1] = -k/m * u[1] - fd_friction/m

    elseif ut[1] < 0.0

        utt[1] = -k/m * u[1] + fd_friction/m

    end

end;

#' Define the normal force of the building.
g = 9.8 #m/s^2
normal_force = m * g 

#' Assume there is a friction force working to dissipate energy in the building (maybe from the shear walls?).  Define this force as a ratio of the total normal force of the building (maybe not so physically meaningful).
fd_friction = 0.1 * normal_force

#' Package up the important physical parameters.
p = [k, m, fd_friction];

#' Build our model. 
prob = SecondOrderODEProblem(equation_of_motion_with_friction_damping, ut_o, u_o, (0.0, 1.0), p);

#' And solve.  
solution = solve(prob, DPRKN6(), tstops=0:0.001:1.0);

#' And plot.
using Plots
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "lateral disp. u [m]")

