using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading
# using Plots

#Calculate the response of a two DOF dynamic system.

#two masses
m1 = 1.0
m2 = 1.0

#two springs
k1 = 1.0
k2 = 1.0

#two viscous dampers
c1 = 0.0
c2 = 0.0

#define mass, stiffness, and damping matrices
M = [m1 0
     0 m2]
K = [k1+k2 -k2
     -k2   k2]
C = [c1+c2 0
     -c2  c2]

#define arbitrary forcing functions
#at DOF 1
t1 = 0:0.01:10 #seconds
a1 = 0.0
ω1 = 1.0
pt1 = a1*sin.(ω1*t1)
pt1Spline = Spline1D(t1,pt1)

#at DOF 2
t2 = 0:0.01:3 #seconds
a2 = 0.0
ω2 = 1.0
pt2 = a2*cos.(ω2*t2)
pt2Spline = Spline1D(t2,pt2)

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt2Spline = p

    ddu[:,1] = -inv(M)*C*du[:,1]-inv(M)*K*u[:,1]

end
 #                                    u_dot0      u_0     trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.], [-1.; 1.], (0.,10.),(M, C, K, pt1Spline, pt2Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u1=(x->x[3]).(sol.u)
u2=(x->x[4]).(sol.u)

t=sol.t


#animate the 2 DOF dynamic system
using Makie


#offset mass1 and mass2 from each other

u1=u1 .+ 3.0
u2=u2 .+ 6.0


#define initial position for mass-spring animation
XYPos1=Node((u1[1],0))
XYPos2=Node((u2[1],0))

#define window and limits
scene1=Scene(resolution=(1000,500))
Xmin=minimum([u1;u2])-1.0
Ymin=-1
Xmax=maximum([u1;u2])+1.0
Ymax=1
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#show mass
scatter!(scene1, lift(xy->Point2f0[xy],XYPos1), show_axis = false, marker = [:rect], limits=limits, color= :red, markersize=0.5)
scatter!(scene1, lift(xy->Point2f0[xy],XYPos2), show_axis = false, marker = [:rect], limits=limits, color= :blue, markersize=0.5)
lines!([0, 0], [-1, 1], color= :black)
lines!([3, 3], [-1, 1], color= :red)
lines!([6, 6], [-1, 1], color= :blue)

#animate
#adjust step in loop to match actual time with animation time
record(scene1, "animation.mp4", range(1, stop = length(u1), step=3)) do i

    # val, looptime, bytes, gctime, memallocs=@timed begin  #measure time for each loop
    XYPos1[]=(u1[i],0)
    XYPos2[]=(u2[i],0)

end
