using OrdinaryDiffEq  #ODE solver
using DifferentialEquations
using Plots
using Dierckx  #for interpolating an aribtrary loading

#write equation of motion function
#        m,g,d,pt1,pt2,pt3,pwt1,pwt2,pwt3 = p
#    ddu[1] = (pt1+pwt1)/m
#    ddu[3] = (pt3+pwt3)/m
#    ddu[4] = (-pt3-pwt3)/m*d
#    ddu[5] = 0
#    ddu[6] = (-pt1-pwt1)/m*d
#end

#write equation of motion function
function drone!(du,u,p,t)
        m,d,pt1Spline,pt2Spline,pt3Spline,pwt1Spline,pwt2Spline,pwt3Spline = p
        pt1=pt1Spline(t[1])*launch(t[1])*land(t[1])
        pt2=pt2Spline(t[1])*land(t[1])-1000*landland(t[1])
        pt3=pt3Spline(t[1])*launch(t[1])*land(t[1]) #engine forces acting in the x,y,z direction
        pwt1=pwt1Spline(t[1])*launch(t[1])*land(t[1])
        pwt2=pwt2Spline(t[1])*launch(t[1])*land(t[1])
        pwt3=pwt3Spline(t[1])*launch(t[1])*land(t[1]) #wind forces acting in the x,y,z direction
    du[1] = u[2]
    du[2] = (pt1+pwt1)/m
    du[3] = u[4]
    du[4] = (pt2+pwt2)/m
    du[5] = u[6]           # u[1],u[3],u[5] represents linear displacement in x, y, z directions
    du[6] = (pt3+pwt3)/m   # u[2],u[4],u[6] represents linear velocity in x, y, z directions
    du[7] = u[8]
    du[8] =  (-pt3-pwt3)/(m*d)
    du[9] = u[10]
    du[10] = 0
    du[11] = u[12]              # u[7],u[9],u[11] represents angular displacement in x, y, z directions
    du[12] = (-pt1-pwt1)/(m*d)  # u[8],u[10],u[12] represents angular velocity in x, y, z directions
end


#estimate parameters
m = 1000 # weight of the drone
d = 1 # distance between the mass point(center of the roller) and the position forces applied to (center of the drone body)

#estimate wind forces
t0 = 0:0.1:100 #seconds
a0 = 1000
ω0 = 1
α = π/2 # angle between the x-axis and the direction of wind force among xy-plane
β = π/2 # angle between the xy-plane and the direction of wind force
pwt = a0*sin.(ω0*t0) # wind force
pwt1Spline = Spline1D(t0,pwt*cos(α)*cos(β)) # x-axis componet
pwt2Spline = Spline1D(t0,pwt*sin(β)) # y-axis componet
pwt3Spline = Spline1D(t0,pwt*sin(α)*cos(β)) # z-axis componet

#define arbitrary forcing functions
#at x direction
t1 = 0:0.1:100 #seconds
a1 = 1000
ω1 = 1
p1 = a1*sin.(ω1*t1) # engine force in x-direction
pt1Spline = Spline1D(t1,p1) # interpolated function of time and p1


#at y direction
function launch(t)
   0.5 * (sign(t-20) + 1) # ensure no force except vertical force(pt2) acts before 20s
end
function land(t)  # ensure no force except vertical force acts after 80s
   0.5 * (sign(80-t)+1)
end
function landland(t) # ensure the vertical force changes to its opposite(-pt2) acts after 83.6s
   0.5 * (sign(t-83.6) + 1)
end
t2 = 0:0.1:100 #seconds
a2 = 1000
ω2 = 1
p2 = a2*sin.(ω2*t2) # engine force in y-direction
pt2Spline = Spline1D(t2,p2)



#at z direction
t3 = 0:0.1:100 #seconds
a3 = 0
ω3 = 0
p3 = a3*sin.(ω3*t3) # engine force in z-direction
pt3Spline = Spline1D(t3,p3)



udot=[0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
tspan=(0.0,100.0)
prob = ODEProblem(drone!,udot,tspan,(m,d,pt1Spline,pt2Spline,pt3Spline,pwt1Spline,pwt2Spline,pwt3Spline))
sol = solve(prob) # solve the ODE function

plot(sol) # plot all of the displacements and velocities
