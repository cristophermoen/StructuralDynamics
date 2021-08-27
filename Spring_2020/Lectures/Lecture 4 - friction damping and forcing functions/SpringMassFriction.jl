#define Julia ODE package that we will use as solver
using OrdinaryDiffEq

#define constants
m=3000000    #kg
g=9.8  #m/sec^2
Tn=0.5  #seconds
fn=1/Tn
ωn=2*pi*fn   #radians/second

k=ωn^2*m    #N/m

ζ=0.05      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

μ=0.1  #friction coefficient
N=m*g   #normal force
F=μ*N   #friction force

#calculate frequency and period
ωn=sqrt(k/m)  #radians per second
fn=ωn/(2*pi)  #cycles per second
Tn=1/fn       #period, seconds

#forcing function
a=0 #N
ω=0  #radians/sec.

#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, a, ω, F  = p

    #harmonic forcing function
    pt=a*sin(ω*t[1])

    #if friction force is bigger than spring force, stop solution
    #friction damping
    if du[1]>0   #velocity
        SignF=1
    else
        SignF=-1
    end

    ddu[1] = -c/m*du[1] - k/m*u[1] + pt/m -F*SignF/m
end

#solver controls
dt=0.001 #seconds
totalt=10 #seconds

#define initial conditions
uo=0.1  #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m, k, c, a, ω, F))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)   #angular velocity
t=sol.t

#plot results
using Plots    #start Julia plot package
plot(t,u,linewidth=1,
    xaxis="t [sec.]",yaxis="u [meters]", legend=false)
