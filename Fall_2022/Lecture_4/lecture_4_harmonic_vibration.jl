# Lecture 4
# EN 560.630 Fall 2022
# Department of Civil and Systems Engineering
# Johns Hopkins University

using DifferentialEquations, Plots

m = 100000.0 #kg
k = 30000000.0 #kN/m

ωn = sqrt(k/m)

ζ = 0.05
c_cr = 2*m*ωn
c = ζ * c_cr


po = 100000.0  #kN
ω = 0.3 * ωn

#package up parameters
p = [m, k, c, po, ω]

# Define an equation of motion.
function equation_of_motion(utt, ut, u, p, t)
	
    m, k, c, po, ω = p

    utt[1] = - c/m * ut[1] -k/m * u[1] - (po/m)* sin(ω*t[1])

end

# Define initial conditions.
ut_0 = [0.0] #m/s 
u_0 = [0.5] #m

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 10.0
dt = 0.01

# Build our model. 
prob = SecondOrderODEProblem(equation_of_motion, u_0, ut_0, (t_start, t_end), p);

#And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#And plot.
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [m]")

ut = (x->x[1]).(solution.u)  #velocity
plot(solution.t, ut, legend = false, xlabel="time [seconds]", ylabel = "ut [mm/sec]")




