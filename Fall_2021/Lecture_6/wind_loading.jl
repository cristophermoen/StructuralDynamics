#' # Lecture 6 - Simulating a wind loading 

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Wind pressure varies over time, both in magnitude and in its frequency.

#' Let's model a random wind pressure. 

#' Define a mean wind velocity.
μ_wind_velocity = 5.0  #m/s

#' Define the standard deviation of the wind velocity.
σ_wind_velocity = 5.0 #m/s

#' Calculate wind pressure from wind velocity where velocity is in m/s.
function calculate_wind_pressure_in_metric(V)

    #Convert m/s to mph.
    V = V / 0.44704

    pressure = 0.00256 * V^2

    #Convert pressure from psf to N/m^2.
    pressure = pressure * 47.88

    return pressure

end

μ_wind_pressure = calculate_wind_pressure_in_metric(μ_wind_velocity)
σ_wind_pressure = calculate_wind_pressure_in_metric(σ_wind_velocity)

#' Define a wind area over which the pressure acts.
A_wind = 1.0  #m^2 

#' Calculate the wind force statistics.  
μ_wind_force = μ_wind_pressure * A_wind
σ_wind_force = σ_wind_pressure * A_wind

#' Define a Gaussian distribution for the wind force.
using Distributions
wind_pdf = Normal(μ_wind_pressure, σ_wind_pressure)

#' Define an equation of motion including a random wind force.
function equation_of_motion(utt, ut, u, p, t)
	
    k, m, c, wind_load_function = p

    utt[1] = -k/m * u[1] - c/m * ut[1] - wind_load_function(t[1])/m

end;

#' Define the mass of the tower.
m = 1000 #kg

#' Define the stiffness of the tower
k = 100  #N/m

#' Define viscous damping parameters.
ωn = sqrt(k/m)   #radians/sec
c_cr = 2 * ωn * m  #kg/sec
ξ = 0.05
c = ξ * c_cr  #kg/sec

#' We don't need an initial velocity or initial displacement to get the tower going since we are using a forcing function.
ut_o = [0.0] 
u_o = [0.0]

#' Define the time range over which the model should run.
t_start = 0.0
t_end = 10.0
dt = 0.001

#' Calculate the wind load distribution.
t_dist_start = 0.0 #sec
dt_dist = 0.1 #sec
t_dist_end = 10.0 #sec
wind_load_over_time = rand(wind_pdf, length(t_dist_start:dt_dist:t_dist_end))

#' Define this p(t) with a continous function.
using Dierckx
wind_load_function = Spline1D(t_dist_start:dt_dist:t_dist_end, wind_load_over_time)

#' Package up the important physical parameters.
p = [k, m, c, wind_load_function];

#' Build our model. 
prob = SecondOrderODEProblem(equation_of_motion, ut_o, u_o, (t_start, t_end), p);

#' And solve.  
solution = solve(prob, DPRKN6(), tstops=t_start:dt:t_end);

#' And plot.
using Plots
u = (x->x[2]).(solution.u)  #displacement
plot(solution.t, u, legend = false, xlabel="time [seconds]", ylabel = "u [m]")




