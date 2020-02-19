#define Julia ODE package that we will use as solver
using OrdinaryDiffEq

#define constants
m=3000000    #kg
g=9.8  #m/sec^2
Tn=0.5  #seconds
fn=1/Tn
ωn=2*pi*fn   #radians/second

k=ωn^2*m    #N/m

ζ=0.00      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

μ=0.10  #friction coefficient
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

    #friction damping
    if du[1]>0
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

#animate u vs. friction force plot with a red point follower

using Makie

#define friction force at every time t in the simulation
Ft=F*sign.(du)

#tuple of x-y follower position, the Node command automatically updates later
#in the record function
XYPos=Node((u[1],Ft[1]))

#define window and limits
scene=Scene(resolution=(1000,500))
Xmin=minimum(u)
Ymin=minimum(Ft)
Xmax=maximum(u)
Ymax=maximum(Ft)
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#plot the whole response
lines!(scene, u,Ft, limits=limits, color= :blue)

#show pointer
scatter!(scene, lift(xy->Point2f0[xy],XYPos), marker = [:circle], limits=limits, color= :red, markersize=100000)

#animate
record(scene, "uFtAnimation.mp4", range(1, stop = length(u), step=1)) do i
     XYPos[]=(u[i],Ft[i])
end
