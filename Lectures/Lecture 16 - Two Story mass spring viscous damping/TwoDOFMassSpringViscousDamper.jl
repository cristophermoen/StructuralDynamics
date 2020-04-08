using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading


#Calculate the free vibration response of a two story shear building

#floors
m1 = 8000/2.281  #kg   #first floor, 20 psf dead load
m2 = 4000/2.281  #kg   #second floor (or roof), 10 psf dead load

#columns
E=2E+11  #N/m^2
Ic=723/12^4/3.281^4  #m^4   #W14x68 steel column, strong axis
L=13/3.281  #m  column height
k1 = 12*E*4*Ic/L^3
k2 = 12*E*2*Ic/L^3

#two viscous dampers
c1 = 0.0   #units are kg/sec
c2 = 0.0

#define mass, stiffness, and damping matrices
M = [m1 0
     0 m2]
K = [k1+k2 -k2
     -k2   k2]
C = [c1+c2 -c2
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

    pt1=pt1Spline(t[1])
    pt2=pt2Spline(t[1])
    pt=[pt1; pt2]

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0      u_0     trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.], [1.0; 1.0], (0.,10.),(M, C, K, pt1Spline, pt2Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u1=(x->x[3]).(sol.u)
u2=(x->x[4]).(sol.u)

t=sol.t


#animate the two story shear building
#more realistic now, with columns and beams

using Makie

#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum([u1;u2])-4.0
Ymin=0-4.0
Xmax=maximum([u1;u2])+4.0
Ymax=2*L+4.0
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)




#define initial position
XYPos1=Node(u1[1])
XYPos2=Node(u2[1])


#draw floors
#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L+L/16], [-L/2+xy, L+L/16], [-L/2+xy, L-L/16], [L/2+xy, L-L/16]], XYPos1), color = :lightgray, show_axis = false, limits=limits)
#second floor
poly!(scene1, lift((xy1,xy2)->Point2f0[[L/2+xy2, 2*L+L/16], [-L/2+xy2, 2*L+L/16], [-L/2+xy2, 2*L-L/16], [L/2+xy2, 2*L-L/16]], XYPos1, XYPos2), color = :lightgray, show_axis = false, limits=limits)

#draw ground and building centerline
lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)

#draw masses

scatter!(scene1, lift(xy->Point2f0[(xy,L)],XYPos1), marker = [:circle], limits=limits, color= :red, markersize=1)  #first floor mass
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
