
#' # Homework 2 solution - Steel beam floor vibration from party people jumping, with Lazan material damping 

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' There is a party, and people are jumping up and down to the beat of [House of Pain "Jump Around"](https://www.youtube.com/watch?v=XhzpxjuwZy0).   The party floor is vibrating.   Let's model it.

#' This is the floor.  It is supported by steel W-shape beams.

#' ![Figure 1](/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Homework_solutions/HW1/coord_sys_schematic.png)

#' The plan is to model this dynamic system as two degrees of freedom, one DOF for the beam, and one DOF for the people jumping up and down.

#' Units are kg, meters, seconds, and Newtons.

#' Define the beam span.
L = 20 #m

#' Define the beam spacing.
beam_spacing = 2.0 #m

#' Define the weight of the people on the floor.  Start with an assumed pressure of 100 psf which comes from [ASCE 7-16](https://www.asce.org/publications-and-news/asce-7).   

#' Convert this to metric.
people_pressure = 100 * 47.88  #N/m^2

#' Convert to force over a tributary area on a beam.
people_force = people_pressure * L * beam_spacing  #N

#' Define the acceleration of gravity.
g = 9.81 #m/s^2

#' Convert force to a mass.
people_mass =  people_force/g  #kg

#' Define the steel beam cross-section.  Use a [W14x22](https://engineering-database.com/structural-members/beams/wide-flange-beams/W14x22-a572-50) 
I_beam = 199 #in^4
A_beam = 6.49 #in^4

#' Convert to metric.
I_beam = I_beam * 0.0254^4  #m^4
A_beam = A_beam * 0.0254^2  #m^2

#' Define the density of steel.
γ_steel = 8050  #kg/m^3

#' Calculate the mass of the beam.
beam_mass = A_beam * γ_steel

#' Assume the floor deck that the beam is supporting weights 20 psf.
floor_deck_weight = 20 #psf

#' Convert to metric.
floor_deck_weight = floor_deck_weight * 47.88 #psf

#' Convert to force over a tributary area on a beam.
floor_deck_force = floor_deck_weight * L * beam_spacing  #N

#' Convert force to a mass.
floor_deck_mass =  floor_deck_force/g  #kg

#' Define the total mass of the floor.
floor_mass = beam_mass + floor_deck_mass

#' Approximate the floor stiffness by assuming the beam is loaded uniformly.  The midspan beam deflection is $\Delta = 5wL^4/(384EI)$ and the lumped beam force is $P = wL$.   Use then $k=P/\Delta$ with a unit distributed load.

#' Define the elastic modulus for steel.
E_steel = 1.99948e+11 #N/m^2

Δ = 5*1.0*L^4/(384 * E_steel * I_beam)
P = 1.0 * L
k_beam = P / Δ 

#' Now let's work on the equations of motion.

#' There are two degrees of freedom now, so mass and stiffness will be represented by matrices.

#' The simulation starts when the song starts.  At that moment everyone is standing on the floor, and then they start jumping.

M = [people_mass 0.0
    0.0 floor_mass + people_mass]

#' Use a contact model for the people jumping.  When the displacement u_people becomes equal to u_beam then k_contact is a big number (simulating contact).   
k_contact = 100 * k_beam

K = [k_contact 0.0
     0.0  k_beam + k_contact]    

P = [people_mass * g
     floor_mass * g]  


p = [M, K, P]

function equation_of_motion(utt, ut, u, p, t)

    M, K, P = p

    utt[:,1] = +inv(M)*K*u[:,1] + inv(M)*P

end
using DifferentialEquations
 #                                                  u_dot0      u_0     trange
prob = SecondOrderODEProblem(equation_of_motion, [0.; 0.], [0.; 0.], (0.,1.),(M, K, P))
sol = solve(prob, DPRKN6(),tstops=0.:0.01:1.)


u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u1=(x->x[3]).(sol.u)
u2=(x->x[4]).(sol.u)

t=sol.t



using Plots

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

n = 150
t = range(0, 2π, length = n)
x = sin.(t)
y = cos.(t)

anim = @animate for i ∈ 1:n
    circleplot(x, y, i)
end
gif(anim, "anim_fps15.gif", fps = 15)


using CairoMakie, Makie

x = range(0, 10, length=100)
y = sin.(x)
lines(x, y)


# using GLMakie
# using Makie
# using GeometryBasics


# time = Node(0.0)

# xs = range(0, 7, length=40)

# ys_1 = @lift(sin.(xs .- $time))
# ys_2 = @lift(cos.(xs .- $time) .+ 3)

# fig = lines(xs, ys_1, color = :blue, linewidth = 4,
#     axis = (title = @lift("t = $(round($time, digits = 1))"),))
# scatter!(xs, ys_2, color = :red, markersize = 15)

# framerate = 30
# timestamps = range(0, 2, step=1/framerate)

# record(fig, "time_animation.mp4", timestamps;
#         framerate = framerate) do t
#     time[] = t
# end




#define window and limits
# scene1=Scene(resolution=(1000,1000))
# Xmin=minimum([u1;u2])-4.0
# Ymin=0-4.0
# Xmax=maximum([u1;u2])+4.0
# Ymax=2*L+4.0
# Xrange=Xmax-Xmin
# Yrange=Ymax-Ymin
# limits = FRect(Xmin, Ymin, Xrange, Yrange)




#define initial position
XYPos1=Node(u1[1])
XYPos2=Node(u2[1])
  
  
#   #draw floors
#   #first floor
#   poly!(scene1, lift(xy->Point2f0[[L/2+xy, L+L/16], [-L/2+xy, L+L/16], [-L/2+xy, L-L/16], [L/2+xy, L-L/16]], XYPos1), color = :lightgray, show_axis = false, limits=limits)
#   #second floor
#   poly!(scene1, lift((xy1,xy2)->Point2f0[[L/2+xy2, 2*L+L/16], [-L/2+xy2, 2*L+L/16], [-L/2+xy2, 2*L-L/16], [L/2+xy2, 2*L-L/16]], XYPos1, XYPos2), color = :lightgray, show_axis = false, limits=limits)
  
#   #draw ground and building centerline
#   lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
#   lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)
  
  #draw masses
  
  fig, ax = scatter(lift(xy->Point2f0[(xy,L)],XYPos1), marker = [:circle], limits=(0, 30, 0, 30), color= :red, markersize=1)  #first floor mass
  scatter!(scene1, lift((xy1,xy2)->Point2f0[(xy2,2*L)],XYPos1, XYPos2), marker = [:circle], limits=limits, color= :red, markersize=0.5)  #second floor mass
  
  
  #draw columns
  #use shape function and coordinate system from
  #https://www.youtube.com/watch?v=K7lfyAldj5k
  
  function BeamShape(q1,q2,q3,q4,L,x, offset)
  
      a0=q1
      a1=q2
      a2=1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
      a3=1/L^3*(2*q1+q2*L-2*q3+q4*L)
  
      w=a0 .+a1.*x .+a2.*x.^2 .+a3.*x.^3 .+offset
  
  end
  
  x=0:L/10:L
  q1=0
  q2=0
  q4=0
  
  lines!(lift(xy->BeamShape(q1,q2,xy,q4,L,x, L/2),XYPos1),x, linewidth=12)   #first level right column
  lines!(lift(xy->BeamShape(q1,q2,xy,q4,L,x, -L/2),XYPos1),x, linewidth=12) #first level left column
  
  lines!(lift((xy1, xy2)->BeamShape(q1,q2,-xy1+xy2,q4,L,x,L/2+xy1),XYPos1,XYPos2),x .+L, linewidth=6) #second level right column
  lines!(lift((xy1, xy2)->BeamShape(q1,q2,-xy1+xy2,q4,L,x,-L/2+xy1),XYPos1,XYPos2),x .+L, linewidth=6) #second level left column
  
  #animate
  
  record(scene1, "animation.mp4", range(1, stop = length(u1), step=1)) do i
  
      if i==1
          sleep(1)
      end
  
      XYPos1[]=(u1[i])
      XYPos2[]=(u2[i])
      sleep(1/24)
  
  end

