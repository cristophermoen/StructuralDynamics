using DifferentialEquations
using Plots
using Dierckx  #for interpolating an aribtrary loading

#write equation of motion function
function drone!(du,u,p,t)
        m,Imx,Imz,d,pt1Spline,pt2Spline,pt3Spline,pwt1Spline,pwt2Spline,pwt3Spline = p
        pt1=pt1Spline(t[1])
        pt2=pt2Spline(t[1])
        pt3=pt3Spline(t[1])
        pwt1=pwt1Spline(t[1])
        pwt2=pwt2Spline(t[1])
        pwt3=pwt3Spline(t[1])
    du[1] = u[2]
    du[2] = (pt1+pwt1)/m
    du[3] = u[4]
    du[4] = (pt2+pwt2)/m
    du[5] = u[6]
    du[6] = (pt3+pwt3)/m
    du[7] = u[8]
    du[8] =  (-pt3-pwt3)/Imx*d
    du[9] = u[10]
    du[10] = 0
    du[11] = u[12]
    du[12] = (-pt1-pwt1)/Imz*d
end


#estimate parameters
m = 1000
g = 9.8
d = 1
c = 1

#estimate wind forces
t0 = 0:0.1:100 #seconds
a0 = 1
ω0 = 1
α = π/2
β = π/2
pwt = a0*sin.(ω0*t0)
pwt1Spline = Spline1D(t0,pwt*cos(α)*cos(β))
pwt2Spline = Spline1D(t0,pwt*sin(β))
pwt3Spline = Spline1D(t0,pwt*sin(α)*cos(β))

#define arbitrary forcing functions
#at x direction
t1 = 0:0.1:100 #seconds
a1 = 1
ω1 = 1
p1 = a1*sin.(ω1*t1)
pt1Spline = Spline1D(t1,p1)


#at y direction
t2 = 0:0.1:100 #seconds
a2 = 1
ω2 = 1
p2 = a2*sin.(ω2*t2)
pt2Spline = Spline1D(t2,p2)


#at z direction
t3 = 0:0.1:100 #seconds
a3 = 0
ω3 = 0
p3 = a3*sin.(ω3*t3)
pt3Spline = Spline1D(t3,p3)



udot=[0.1,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.0,0.0,0.0,0.0]
tspan=(0.0,100.0)
prob = ODEProblem(drone!,udot,tspan,(m,Imx,Imz,d,pt1Spline,pt2Spline,pt3Spline,pwt1Spline,pwt2Spline,pwt3Spline))
sol = solve(prob, tstops=0:0.1:100)


t=sol.t
ux_dot=(x->x[1]).(sol.u)
ux=(x->x[2]).(sol.u)
uy_dot=(x->x[3]).(sol.u)
uy=(x->x[4]).(sol.u)
uz_dot=(x->x[5]).(sol.u)
uz=(x->x[6]).(sol.u)
θx_dot=(x->x[7]).(sol.u)
θx=(x->x[8]).(sol.u)
θy_dot=(x->x[9]).(sol.u)
θy=(x->x[10]).(sol.u)
θz_dot=(x->x[11]).(sol.u)
θz=(x->x[12]).(sol.u)



plot(sol)
