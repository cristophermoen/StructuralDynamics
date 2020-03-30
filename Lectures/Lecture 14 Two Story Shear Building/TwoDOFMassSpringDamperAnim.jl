using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading


#Calculate the free vibration response of a two story shear building

#floors
m1 = 8000/2.281  #kg   #first floor, 20 psf dead load over a 20x20 ft floorplan
m2 = 4000/2.281  #kg   #second floor (or roof), 10 psf dead load

#columns
E=2E+11  #N/m^2
Ic=723/12^4/3.281^4  #m^4   #W14x68 steel column, strong axis
L=13/3.281  #m  column height
k1 = 12*E*4*Ic/L^3
k2 = 12*E*2*Ic/L^3

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
  prob = SecondOrderODEProblem(mdof, [0.; 0.], [1.0; 0.], (0.,10.),(M, C, K, pt1Spline, pt2Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u1=(x->x[3]).(sol.u)
u2=(x->x[4]).(sol.u)

t=sol.t


#animate the two story building

using Makie

#define initial position
XYPos1=Node((u1[1],L))
XYPos2=Node((u2[1],2*L))

#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum([u1;u2])-4.0
Ymin=0-4.0
Xmax=maximum([u1;u2])+4.0
Ymax=2*L+4.0
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#show mass
scatter!(scene1, lift(xy->Point2f0[xy],XYPos1), show_axis = false, marker = [:rect], limits=limits, color= :red, markersize=1.0)
scatter!(scene1, lift(xy->Point2f0[xy],XYPos2), show_axis = false, marker = [:rect], limits=limits, color= :blue, markersize=0.5)
lines!([Xmin, Xmax], [0, 0], color= :black)
lines!([0, 0], [Ymin, Ymax], color= :black)

#animate
#adjust step in loop to match actual time with animation time
record(scene1, "animation.mp4", range(1, stop = length(u1), step=10)) do i

    # val, looptime, bytes, gctime, memallocs=@timed begin  #measure time for each loop
    XYPos1[]=(u1[i],L)
    XYPos2[]=(u2[i],2*L)

end
