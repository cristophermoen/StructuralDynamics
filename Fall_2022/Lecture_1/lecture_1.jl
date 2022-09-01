# Lecture 1 - solving equation of motion
# EN 560.630 Fall 2022
# Department of Civil and Systems Engineering
# Johns Hopkins University

using DifferentialEquations, Plots

# Define an equation of motion.
function equation_of_motion(utt, ut, u, p, t)
	
    m = p

    utt = 0.0

end

# Define mass.
m = 1000 #kg

# Define initial conditions.
ut_o = [1.0] #m/s 
u_o = [0.0] #m

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 10.0
dt = 0.01

# Build our model. 
prob = SecondOrderODEProblem(equation_of_motion, ut_o, u_o, (t_start, t_end), (m));

#And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#And plot.
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [m]")

ut = (x->x[1]).(solution.u)  #velocity
plot(solution.t, ut, legend = false, xlabel="time [seconds]", ylabel = "u [m]")




