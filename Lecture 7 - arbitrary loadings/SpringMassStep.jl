#define Julia ODE package that we will use as solver
using OrdinaryDiffEq
using Dierckx  #this is for interpolating an aribtrary loading

#define constants
m=300    #kg
g=9.8  #m/sec^2
Tn=0.5  #seconds
fn=1/Tn
ωn=2*pi*fn   #radians/second

k=ωn^2*m    #N/m

ζ=0.05      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

μ=0.00  #friction coefficient
N=m*g   #normal force
F=μ*N   #friction force

#calculate frequency and period
ωn=sqrt(k/m)  #radians per second
fn=ωn/(2*pi)  #cycles per second
Tn=1/fn       #period, seconds

#forcing function
a=0 #N
ω=0  #radians/sec.

#define arbitrary function
#in this case a rectangular step function

function RectangularStepLoading(StepMagnitude,StepTimeStart,StepTimeEnd,t)

    if (t>=StepTimeStart) & (t<=StepTimeEnd)
        fStep=StepMagnitude
    else
        fStep=0.0
    end

    return fStep
end

StepMagnitude=100.0
StepTimeStart=0.0
StepTimeEnd=4.0
tArb=0:0.001:10
fArb=RectangularStepLoading.(StepMagnitude,StepTimeStart,StepTimeEnd,tArb)

plot(tArb,fArb)

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)

#check that spline is matching curve
fArbCheck=spl.(tArb)
plot!(tArb,fArbCheck)

#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, a, ω, F, spl  = p

    #harmonic forcing function
    ptForcing=a*sin(ω*t[1])

    #interpolate abitrary force
    ptArb = spl(t[1])

    #friction damping
    if du[1]>0
        SignF=1
    else
        SignF=-1
    end

    ddu[1] = -c/m*du[1] - k/m*u[1] + ptForcing/m -F*SignF/m + ptArb/m
end

#solver controls
dt=0.0001 #seconds
totalt=10 #seconds

#define initial conditions
uo=0.0  #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m, k, c, a, ω, F, spl))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)   #angular velocity
t=sol.t

#plot results
using Plots    #start Julia plot package
plot(t,u,linewidth=1,
    xaxis="t [sec.]",yaxis="u [meters]", legend=false)
