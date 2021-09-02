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

#' Define the tennis ball spring constant.  Assume the ball has an elastic stiffness of 10 lbs / 0.7 in.   
k_ball = 0.01* 4.448 * 1000 / (0.7 * 25.4/1000); #N/m

#' When the ball is not in contact with the door, set the spring constant to zero.
k_no_contact = 0.0;

#' The equation of motion when the ball is in contact with the door is $m \ddot u + k(u - u_{impact}) = 0$.   The spring acts in a coordinate system relative to the center of the tennis ball, and so the global ball movement $u$ needs to be subtracted by the distance to the door.  Let's call this term $u_{impact}$.'
u_impact = (x_door - d_ball/2);


#' Let's start defining our model in the Julia language.

#' Initialize the system stiffness.
k = k_no_contact;

#' Package up the important physical parameters.
p = [k, m, u_impact];

#' Define the equation of motion.  The acceleration is $u_{tt}$, the velocity is $u_{t}$, the displacement is $u$, all at time $t$.  When $k = 0$, the equation of motion defaults to $m \ddot u = 0$ which is the case when the ball is not in contact with the door.
function equation_of_motion(utt, ut, u, p, t)
	
    k, m, u_impact = p

	utt[1] = -k/m * (u[1] .- u_impact)

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

#' And solve.
solution = solve(prob, DPRKN6(), callback=cbs);

#' And plot.  
using Plots
plot(solution)

#' This model obviously doesn't work well.   Let's try in a different script to add some energy dissipation when the ball hits the door.






