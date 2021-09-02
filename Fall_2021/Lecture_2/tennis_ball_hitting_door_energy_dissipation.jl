#' # Lecture 2 - modeling a tennis ball thrown at the classroom door 

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' Use Julia language package DifferentialEquations.jl to simulate the motion of a tennis ball bouncing off our classroom door.
using DifferentialEquations

#' Units are kg, meters, seconds, and Newtons.

#' Assume a unidirectional coordinate system in $x$.  

#' Define the distance from the person throwing the ball to the door.
x_door = 1.0 #m;

#' Define the tennis ball diameter.
d_ball = 0.068 #m;

#' Assume the tennis ball begins traveling with an initial velocity imparted by the thrower.
ut_o = [1.0];  #m/s

#' The tennis ball starts at zero displacement.
u_o = [0.0];

#' The equation of motion when the ball is released is $m \ddot u = 0$ which assumes the ball is traveling at the constant initial velocity $\dot u_o$.

#' Define the mass of the tennis ball.
m = 0.0577; #kg

#' When the tennis ball hits the door, assume the tennis ball deforms like a spring.

#' Define the tennis ball spring constant. Assume the ball has an elastic stiffness of 10 lbs / 0.7 in.   
k_ball = 0.01* 4.448 *1000 / (0.7 * 25.4/1000); #N/m

#' When the ball is not in contact with the door, set the spring constant to zero.
k_no_contact = 0.0;

#' Consider some energy dissipation when the tennis ball hits the door.  This model will assume that some percentage of the ball's kinetic energy is dissipated while in contact with the door.

#' Calculate the ball kinetic energy before it hits the wall.
kinetic_energy = 1/2 * m * ut_o[1]^2

#' Make an assumption for how much the ball deforms when it hits the door.  
max_ball_deformation = 0.005  #m

#' Assume displacement proportional energy dissipation as the ball hits the door. 

#' Define the amount of energy to be dissipated.
E_diss = 0.10 * kinetic_energy  #N-m

#' We will consider this energy dissipation as a force multiplied by the ball's deformation.  This is displacement proportional damping.

#' Half of the energy is dissipated from the time when the door contacts the wall to when the ball's velocity is zero.  When the ball rebounds, assume the other half of the E_diss is consumed.  This is why there is a divide by 2 below.
Fd = E_diss/(1/2 * max_ball_deformation) / 2   #N

#' We need the slope κ of the damping force vs. ball deformation curve for displacement proportional damping.   
κ = Fd/max_ball_deformation  #N/m

#' The equation of motion when the ball is in contact with the door is then $m \ddot u + ku{ball_def} + κ * u_{ball_def} = 0$ when $\dot u > 0$ and $m \ddot u + ku{ball_def} - κ * u_{ball_def} = 0$ when $\dot u < 0$.  This means that for either direction that the ball is moving (towards the door or away from the door) that energy is being dissipated.  You can understand this better by drawing the free body diagram of the ball at different times.  See Chopra Figure 1.5.1 for an example.'

#' Let's start defining our model in the Julia language.

#' Initialize the system stiffness.
k = k_no_contact;

#' Package up the important physical parameters.
p = [k, m, κ, x_door, d_ball];

#' Define the equation of motion.  The acceleration is $u_{tt}$, the velocity is $u_{t}$, the displacement is $u$, all at time $t$.'

function equation_of_motion(utt, ut, u, p, t)
	
    k, m, κ, x_door, d_ball = p

    if k == 0.0

        utt[1] = 0.0   #tennis ball is not in contact with the door

    elseif k > 0.0

        ball_deformation = u[1] .- (x_door - d_ball/2)

        if (ut[1] > 0.0)

    	    utt[1] = -k/m * ball_deformation - κ/m * ball_deformation   #tennis ball is moving towards the door surface

        else

            utt[1] = -k/m * ball_deformation + κ/m * ball_deformation  #tennis ball is moving away from the door surface

        end

    end

end;

#' We need tell our model when $k$ changes.  The first condition is when the ball hits the door, and this occurs when the distance the ball has traveled is greater than or equal to $u_{impact}$.
function ball_hits_door(u,t,integrator) 

    u[2] >= u_impact

end;

#' The action caused by this condition, when the door becomes in contact with the door, is that $k=k_{ball}.  The Julia language sees $k$ as the first value in the $p$ vector.    
function ball_hits_door_action!(integrator)

    integrator.p[1] = k_ball  
   
end;

#' The ball rebounds from the door when the ball location is less than $u_{impact}$.
function ball_leaves_door(u,t,integrator) 

    u[2] <= u_impact

end;

#' And the action is that k goes back to zero.
function ball_leaves_door_action!(integrator)

    integrator.p[1] = 0.0  
   
end;

#' These actions are referred to as _callbacks_ in the Julia language.
cb1 = DiscreteCallback(ball_hits_door, ball_hits_door_action!);
cb2 = DiscreteCallback(ball_leaves_door, ball_leaves_door_action!);

#' When you have multiple actions, you can string them together as a set.
cbs = CallbackSet(cb1, cb2);

#' Now we can build our model.   
prob = SecondOrderODEProblem(equation_of_motion, ut_o, u_o, (0.0, 2.0), p);

#' And solve.   I had to tell the solver to use small timesteps to get an accurate solution, see the tstops command below.
solution = solve(prob, DPRKN6(), callback=cbs, tstops=0:0.001:2.0);

#' And plot.  
using Plots
plot(solution)

#' Calculate the maximum ball deformation.

#' First pull out the displacement vector from the solution object.
u = (x->x[2]).(solution.u)

#' Calculate the maximum displacement.
u_max = maximum(u)

#' This is the maximum ball displacement.  It is consistent with the max_ball_deformation magnitude we assumed when defining the energy dissipation model. 
max_ball_deformation = u_max - u_impact

