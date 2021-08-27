using OrdinaryDiffEq  #define Julia ODE package that we will use as solver
using Dierckx  #this is for interpolating an aribtrary loading
using Plots    #start Julia plot package

#define constants
mtrampoline=5.0 #kg
mperson=60.0 #kg
g=9.8  #m/sec^2

L=5 #width of trampoline, meters

# stiffness k=W/Δ
# base stiffness, assumed from experience
kinitial=3920 #N/M   #EA/L   #this is 60kg person at standstill, trampoline deflects 0.15m
EA=kinitial*L/2


m=mtrampoline+mperson
ωn=sqrt(kinitial/m)
# fn=ωn/(2*π)
# Tn=1/fn

ζ=0.04 #viscous damping ratio. assumed
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

#define a rectangular step function
#assume the jumps steadily and stays on the trampoline for 0.25s
#assume the person jumps 1m into the air. Using kinematic this means the
    #person is in the air for 0.9s. Round this to 1s for simplicity

function RectangularStepLoading(StepMagnitude,time)
    if (time>=0.0) & (time<=0.25)
        fStep=StepMagnitude
    elseif (time>=1.25) & (time<=1.5)
        fStep=StepMagnitude
    elseif (time>=2.5) & (time<=2.75)
        fStep=StepMagnitude
    elseif (time>=3.75) & (time<=4.0)
        fStep=StepMagnitude
    elseif (time>=5.0) & (time<=5.25)
        fStep=StepMagnitude
    elseif (time>=6.25) & (time<=6.5)
        fStep=StepMagnitude
    elseif (time>=7.5) & (time<=7.75)
        fStep=StepMagnitude
    elseif (time>=8.75) & (time<=9.0)
        fStep=StepMagnitude
    elseif (time>=10.0) & (time<=10.25)
        fStep=StepMagnitude
    else
        fStep=0.0
    end
    return fStep
end


StepMagnitude=-600 #N or about 60 kg
tArb=0:0.001:10
fArb=RectangularStepLoading.(StepMagnitude,tArb)

plot(tArb,fArb,title="Rectangular Step Function", legend=false,
    xaxis="t (sec)", yaxis="Step Magnitude (Newton)", color=:red)

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)


#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    mperson, mtrampoline, EA, L, c, spl  = p

    #interpolate abitrary force
    ptArb = spl(t[1])

    if ptArb<=0.0 #if person is on the trampoline
        m=mperson+mtrampoline
    else #if person is in the air
        m=mtrampoline
    end

    #k constantly changes depending on displacement for a cable
    k=8EA/L^3*u[1]^2 #derived from https://www.eng-tips.com/viewthread.cfm?qid=388481

    ddu[1] = -c/m*du[1] - k/m*u[1] + ptArb/m
end

#solver controls
dt=0.01 #seconds
totalt=10 #seconds

#define initial conditions
uo=0.0  #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (mperson, mtrampoline, EA, L, c, spl))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)
t=sol.t

#plot results
plot(t,u,linewidth=1,xaxis="t (sec)",yaxis="u (meters)",title="Trampoline displacement vs time", label="trampoline displacement")

plot!(tArb,fArb/600, label="step function/600",legend=:bottomright,
 title="xi=$ζ, kinitial=$kinitial N/m",color=:red )  #plot step function over displacements to see when person lands
