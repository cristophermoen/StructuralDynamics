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

#Define location of wall.
u_wall = 75.0 #mm

#Define the stiffness of the wall.  Very hard!
k_wall = 10000.0 #N/mm

#Define impact conditions callbacks.

function car_hits_wall(u, t, integrator) 

    u[2] > u_wall

end;

#' The action caused by this condition, when the car hits the wall, is that $k=k_{ball}$.  The Julia language sees $k$ as the third value in the $p$ vector.    
function car_hits_wall_action!(integrator)

    integrator.p[3] = k_wall 
   
end;

#' Car no longer in contact with the wall.
function car_leaves_wall(u, t, integrator) 

    u[2] < u_wall

end;

#' And the action is that $k$ goes back to zero.
function car_leaves_wall_action!(integrator)

    integrator.p[3] = 0.0  
    
end;

#' These actions are referred to as _callbacks_ in the Julia language.
cb1 = DiscreteCallback(car_hits_wall, car_hits_wall_action!);
cb2 = DiscreteCallback(car_leaves_wall, car_leaves_wall_action!);

#' When you have multiple actions, you can string them together as a set.
cbs = CallbackSet(cb1, cb2);


# Initialize k to free driving.
k=0.0


#package up parameters
p = [m, fD, k]

# Define an equation of motion.
function equation_of_motion(utt, ut, u, p, t)
	
    m, fD, k = p

    if k==0.0

        if ut[1] > 0.0     #friction force changes direction with velocity
            utt[1] = -fD/m
        else
            utt[1] = fD/m
        end

    else  #when car is in contact with the wall

        utt[1] = -k/m

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
solution = solve(prob, DPRKN6(), callback=cbs, tstops=t_start:dt:t_end);

#And plot.
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [mm]")

ut = (x->x[1]).(solution.u)  #velocity
plot(solution.t, ut, legend = false, xlabel="time [seconds]", ylabel = "ut [mm/sec]")



