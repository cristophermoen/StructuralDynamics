
#' # Homework 1 solution - Tennis ball dropped from a height, bouncing on the floor

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' Model a tennis ball bouncing on the floor after being dropped at a certain height.

using DifferentialEquations

#' Units are kg, meters, seconds, and Newtons.

#' Assume a unidirectional coordinate system in $y$ where $y=0$ is the ground.

#' Define the tennis ball diameter.
d_ball = 0.068 #m;

#' Assume no initial velocity.
ut_o = [0.0];  #m/s

#' The tennis ball starts at a displacement of 1 m from the ground.
u_o = [2.0];

#' Define gravity acceleration constant
g = -9.8 #m^2/s;

#' Define the mass of the tennis ball.
m = 0.0577; #kg

#' When the tennis ball hits the ground, assume the tennis ball deforms like a spring.

#' Define the tennis ball spring constant. Assume the ball has an elastic stiffness of 10 lbs / 0.7 in.   
k_ball = 0.01* 4.448 *1000 / (0.7 * 25.4/1000); #N/m

#' When the ball is not in contact with the ground, set the spring constant to zero.
k_no_contact = 0.0;

#' Consider some energy dissipation when the tennis ball hits the ground.  This model will assume that some percentage of the ball's kinetic energy is dissipated while in contact with the ground.

#' Calculate the ball energy before it hits the ground.
potential_energy = m * -g * u_o[1]   #need negative sign here!

#' Make an assumption for how much the ball deforms when it hits the ground the first time.  From preliminary model runs, this deformation magnitude seems to make sense.  20 mm is believable for a 2 meter drop.  This is a little less than half it's diameter.      
max_ball_deformation = 0.020  #m

#' Assume that 50% of the energy in the ball is dissipated on the first bounce.
E_diss = 0.5 * potential_energy  #N-m

#' We will consider this energy dissipation as a force multiplied by the ball's deformation.  This is displacement proportional damping.

#' Half of the energy on the first bounce is dissipated from the time when the ball contacts the ground to when the ball's velocity is zero.  When the ball rebounds, assume the other half of the E_diss is consumed.  This is why there is a divide by 2 below.
Fd = E_diss/(1/2 * max_ball_deformation) / 2    #N

#' We need the slope κ of the damping force vs. ball deformation curve for displacement proportional damping.   
κ = Fd/max_ball_deformation  #N/m

#' Use this same slope κ for subsequent bounces.  

#' Let's start defining our model in the Julia language.

#' Initialize the system stiffness.
k = k_no_contact;

#' Package up the important physical parameters.
p = [k, m, g, κ, d_ball];

#' Define the equation of motion.  The acceleration is $u_{tt}$, the velocity is $u_{t}$, the displacement is $u$, all at time $t$.'

function equation_of_motion(utt, ut, u, p, t)
	
    k, m, g, κ, d_ball = p

    if k == 0.0

        utt[1] = g #tennis ball is not in contact with the ground

    elseif k > 0.0

        ball_deformation = d_ball/2 - u[1]

        if ut[1] < 0.0   #velocity is negative as ball moves towards the ground 

    	    utt[1] = +k/m * ball_deformation + κ/m * ball_deformation + g   

        elseif ut[1] >= 0.0  #velocity is positive as the ball leaves the ground

           utt[1] = k/m * ball_deformation - κ/m * ball_deformation - g 

        end

    end

end;


#' We need tell our model when $k$ changes.  The first condition is when the ball hits the ground.

u_impact = d_ball/2

function ball_hits_ground(u,t,integrator) 

    u[2] <= u_impact

end;

#' The action caused by this condition, when the ball becomes in contact with the ground, is that $k=k_{ball}.  The Julia language sees $k$ as the first value in the $p$ vector.    
function ball_hits_ground_action!(integrator)

    integrator.p[1] = k_ball  
   
end;

#' The ball rebounds from the ground.
function ball_leaves_ground(u,t,integrator) 

    u[2] >= u_impact

end;

#' And the action is that k goes back to zero.
function ball_leaves_ground_action!(integrator)

    integrator.p[1] = 0.0  
    
end;

#' These actions are referred to as _callbacks_ in the Julia language.
cb1 = DiscreteCallback(ball_hits_ground, ball_hits_ground_action!);
cb2 = DiscreteCallback(ball_leaves_ground, ball_leaves_ground_action!);

#' When you have multiple actions, you can string them together as a set.
cbs = CallbackSet(cb1, cb2);

#' Now we can build our model.   
prob = SecondOrderODEProblem(equation_of_motion, ut_o, u_o, (0.0, 5.0), p);

#' And solve.   
solution = solve(prob, DPRKN6(), callback=cbs, tstops=0:0.001:5.0);

#' And plot.  
using Plots
u = (x->x[2]).(solution.u)  #displacement
t = solution.t
plot(t, u, legend=false)
ylims!(-1, 10)

