using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading


#Calculate the free vibration response of a two story shear building

#floors
m1 = 8000/2.281  #kg   #first floor, 20 psf dead load
m2 = 4000/2.281  #kg   #second floor (or roof), 10 psf dead load
m3 = 2000/2.281  #kg   #third floor (or roof), 5 psf dead load

#shear walls per story
E=2E+11/8  #N/m^2   #assume concrete
b=3     #m  width
t=0.2   #m  thickness
A=b*t   #top view cross-sectional area of shear wall
ν=0.20
G=E/(2*(1+ν))
h=13/3.281  #m  story height
k1 = G*A/h
k2 = G*A/h
k3 = G*A/h

#two viscous dampers
c1 = 0.0   #units are kg/sec
c2 = 0.0
c3 = 0.0

#define mass, stiffness, and damping matrices
M = [m1 0 0
     0 m2 0
     0 0  m3]

K = [k1+k2 -k2 0
     -k2   k2+k3 -k3
     0     -k3    k3]
C = [c1+c2 -c2 0
     -c2  c2+c3 -c3
     0     -c3   c3]

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

#at DOF 3
t3 = 0:0.01:3 #seconds
a3 = 0.0
ω3 = 1.0
pt3 = a2*cos.(ω3*t3)
pt3Spline = Spline1D(t3,pt3)

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt2Spline, pt3Spline = p

    pt1=pt1Spline(t[1])
    pt2=pt2Spline(t[1])
    pt3=pt3Spline(t[1])
    pt=[pt1; pt2; pt3]

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0      u_0     trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.], [1.0; 0.0; 0.0], (0.,10.),(M, C, K, pt1Spline, pt2Spline, pt3Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u3_dot=(x->x[3]).(sol.u)
u1=(x->x[4]).(sol.u)
u2=(x->x[5]).(sol.u)
u3=(x->x[6]).(sol.u)

t=sol.t


#animate the three story shear building
#show floors and walls





# **********
using Makie

# function PlotShearBuilding(t, u1, u2, u3, L, m)

L=h


#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum([u1;u2;u3])-4.0
Ymin=0-4.0
Xmax=maximum([u1;u2;u3])+4.0
Ymax=3*L+4.0
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#define initial position
XYPos1=Node(u1[1])
XYPos2=Node(u2[1])
XYPos3=Node(u3[1])


#draw walls

#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L], [-L/2+xy, L], [-L/2, 0], [L/2, 0]], XYPos1), color = :lightgray, show_axis = false, limits=limits, overdraw =:true)
#second floor
poly!(scene1, lift((xy1,xy2)->Point2f0[[L/2+xy2, 2*L], [-L/2+xy2, 2*L], [-L/2+xy1, L], [L/2+xy1, L]], XYPos1,XYPos2), color = :lightgray, show_axis = false, limits=limits)
#third floor
poly!(scene1, lift((xy2,xy3)->Point2f0[[L/2+xy3, 3*L], [-L/2+xy3, 3*L], [-L/2+xy2, 2*L], [L/2+xy2, 2*L]], XYPos2, XYPos3), color = :lightgray, show_axis = false, limits=limits)



#draw floors
#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L+L/16], [-L/2+xy, L+L/16], [-L/2+xy, L-L/16], [L/2+xy, L-L/16]], XYPos1), color = :darkgray, show_axis = false, limits=limits)
#second floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, 2*L+L/16], [-L/2+xy, 2*L+L/16], [-L/2+xy, 2*L-L/16], [L/2+xy, 2*L-L/16]], XYPos2), color = :darkgray, show_axis = false, limits=limits)
#third floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, 3*L+L/16], [-L/2+xy, 3*L+L/16], [-L/2+xy, 3*L-L/16], [L/2+xy, 3*L-L/16]], XYPos3), color = :darkgray, show_axis = false, limits=limits)


#draw ground and building centerline
lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)

#draw masses

scatter!(scene1, lift(xy->Point2f0[(xy,L)],XYPos1), marker = [:circle], limits=limits, color= :red, markersize=1)  #first floor mass
scatter!(scene1, lift(xy->Point2f0[(xy,2*L)],XYPos2), marker = [:circle], limits=limits, color= :red, markersize=0.5)  #second floor mass
scatter!(scene1, lift(xy->Point2f0[(xy,3*L)],XYPos3), marker = [:circle], limits=limits, color= :red, markersize=0.25)  #third floor mass

#animate

record(scene1, "animation.mp4", range(1, stop = length(u1), step=1)) do i

    if i==1
        sleep(1)
    end

    XYPos1[]=(u1[i])
    XYPos2[]=(u2[i])
    XYPos3[]=(u3[i])
    sleep(1/24)

end
