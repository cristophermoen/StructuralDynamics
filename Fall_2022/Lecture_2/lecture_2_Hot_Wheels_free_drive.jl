# Lecture 2
# EN 560.630 Fall 2022
# Department of Civil and Systems Engineering
# Johns Hopkins University

using DifferentialEquations, Plots

ut_0 = 50.0 #mm/sec   #made this faster, maybe underestimated in our experiment
m = 36.0 #grams  #from internet

# all the energy
E_k = 1/2 * m * ut_0^2 

#In our experiment, energy was gone when u=180 mm.
fD = E_k/180.0*1.2  #need to amplify a bit here, to match experiment

#package up parameters
p = [m, fD]

# Define an equation of motion.
function equation_of_motion(utt, ut, u, p, t)
	
    m, fD = p

    utt[1] = -fD/m

    if ut[1]<0.0  #motion stops

        ut[1] = 0.0

    end

end

# Define initial conditions.
ut_0 = [ut_0] #mm/s 
u_0 = [0.0] #m

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 8.0
dt = 0.01

# Build our model. 
prob = SecondOrderODEProblem(equation_of_motion, ut_0, u_0, (t_start, t_end), p);

#And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#And plot.
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [mm]")

ut = (x->x[1]).(solution.u)  #velocity
plot(solution.t, ut, legend = false, xlabel="time [seconds]", ylabel = "ut [mm/sec]")


#Not so bad for rough measurements.

