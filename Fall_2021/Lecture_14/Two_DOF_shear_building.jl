#' # Lecture 13 - Two story shear building earthquake elastic time history simulation

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#You all have enough experience now where we can start modeling dynamic systems with multiple degrees of freedom.   

#' Let start with a two story building subjected to the 1979 Imperial Valley earthquake where the ground motion record is obtained from the [Center for Engineering Strong Motion Data](https://www.strongmotioncenter.org/cgi-bin/CESMD/iqr_dist_DM2.pl?IQRID=ImperialValley79&SFlag=0&Flag=3).

#' Read in ground acceleration data, note dims=(rows,columns) in text file that has 711 lines and 8 columns.
using DelimitedFiles
earthquake_data = readdlm("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_12/Imperial_Valley_1979.txt",Float64,dims=(711,8))

#' Rearrange ground accleration data as a vector.
earthquake_data=transpose(earthquake_data)   #transpose matrix
utt_g=reshape(earthquake_data,(5688,1))  #reshape matrix into a vector, note 8*711=5688
utt_g=vec(utt_g)   #tell Julia to make the data a 1D vector

#' Units come in as cm/sec^2 from strongmotion.org, let's change them to m/sec^2.
utt_g=utt_g ./ 100.0

#' Define time range of the earthquake assuming the first acceleration reading occurs at t=0.0 seconds.  The ground acceleration is provided every 0.010 seconds.
t_eq=collect(range(0,length=5688,stop=5687*0.01))   #total time is 5687*0.01 seconds

#' Plot ground acceleration over time.
using Plots
plot(t_eq,utt_g,linewidth=1,title="Imperial Valley 1979 earthquake",
    xaxis="time, t (sec.)",yaxis="ground acc., m/sec^2",legend=false)

#' Convert the earthquake ground motion from a discrete data set to a continous function.  
using Dierckx
earthquake = Spline1D(t_eq,utt_g)

#' Write the equation of motion.
function mdof(utt, ut, u, p, t)

    M, K, C, earthquake = p

    utt_g = earthquake(t[1])

    #utt[:, 1] = inv(M) * K * u[:, 1] + ....
    utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:,1] - [1.0; 1.0] .* utt_g  
    #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end

#' Describe the building. 

m1 = 10000  #kg   
m2 = 10000  #kg   

k1 = 12000000  #N/m
k2 = 12000000  #N/m 

#' Define mass and stiffness matrices.

M = [m1 0
     0 m2]

K = [k1+k2 -k2
     -k2   k2]

#' Calculate the mode shapes.
using LinearAlgebra

#' Solve for eigenvalues of K*ϕn=wn^2*M*ϕn.
ωn_squared=eigvals(K, M)

#' Calculate the modal natural frequencies.
ωn = sqrt.(ωn_squared)

#' Calculate the mode shape natural frequencies in Hz.
fn = ωn ./(2 * π)

#' Calculate the mode shape natural periods in seconds.
Tn = 1 ./ fn

#' Solve for eigenvectors K*ϕn=wn^2*M*ϕn.
mode_shapes=eigvecs(K,M)

#' Scale the modes so that the maximum modal displacement is 1.0.
num_modes = size(mode_shapes)[2]

scaled_mode_shapes = Matrix{Float64}(undef, size(mode_shapes))

for i = 1:num_modes

     scaling_factor = maximum(abs.(mode_shapes[:,i]))
     scaled_mode_shapes[:,i] = mode_shapes[:,i] ./ scaling_factor
     
end

scaled_mode_shapes

#' Add some damping.

#' Define a viscous damping ratio.
ζ=0.05      #viscous damping ratio
ccr=2 .* [m1; m2] .* ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

C = [c[1] 0.0
     0.0  c[2]]


#' Define the problem.
#'                                   u_dot0      u_0     trange
problem = SecondOrderODEProblem(mdof, [0.; 0.], [0.0; 0.], (0.,60.),(M, K, C, earthquake))

#' Solve.
solution = solve(problem, DPRKN8(),tstops=0:0.01:60.0)

#' Get response.
u1_dot=(x->x[1]).(solution.u)
u2_dot=(x->x[2]).(solution.u)
u1=(x->x[3]).(solution.u)
u2=(x->x[4]).(solution.u)
t=solution.t

#' Plot.
using Plots
plot(t, u1)
plot!(t, u2)

#' Work on an animation.

#' Approximate the ground displacement by integrating utt_g twice.
using QuadGK

num_ground_motion_timesteps = length(t_eq)
ut_g = Array{Float64}(undef, num_ground_motion_timesteps)

for i = 1:num_ground_motion_timesteps
     ut_g[i], err = quadgk(t -> earthquake(t), 0.0, t[i], maxevals=10^4)
end

ground_velocity = Spline1D(t_eq, ut_g)

plot(t_eq, ut_g)
u_g = Array{Float64}(undef, num_ground_motion_timesteps)
for i = 1:num_ground_motion_timesteps
     u_g[i], err = quadgk(t -> ground_velocity(t), 0.0, t[i], maxevals=10^4)
end

plot(t_eq, u_g)

