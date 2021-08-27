#mallika akhtar
#enchilada

#define Julia ODE package that we will use as solver
using OrdinaryDiffEq
using Dierckx  #this is for interpolating an aribtrary loading
using Plots    #start Julia plot package

#from springmassstep.jl
#define constants
m_trampoline=11.3398 #kg - mass of the trampoline
m_person=68.0    #kg -mass of 150 lb person on trampoline
g=9.8  #m/sec^2

L=14*0.3048 #m -14ft length of a trampoline
width=1*0.3048 #m - #1ft affected width (1ft)
A= L*width
E= 456720000 #Pascal
EA=E*A

m=m_trampoline+m_person #kg combined mass of person and trampoline
θ=(m/(2*E*A))^(1/3) #radians - theta=(m/2EA)^1/3

Δ=L/2*tan(θ)

k= m/Δ #stiffness

ωn=sqrt(k/m)   #radians/second
fn=ωn/(2*π)
Tn=1/fn


ζ=0.05      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

#define arbitrary function
#in this case a rectangular step function

function RectangularStepLoading(StepMagnitude,time)

    if (time>=0.0) & (time<=0.20)
        fStep=StepMagnitude

    elseif (time>=1.2) & (time<=1.4)
         fStep=StepMagnitude

    elseif (time>=2.4) & (time<=2.6)
         fStep=StepMagnitude

    elseif (time>=3.6) & (time<=3.8)
         fStep=StepMagnitude

    elseif (time>=4.8) & (time<=5.0)
         fStep=StepMagnitude

    elseif (time>=6.0) & (time<=6.2)
         fStep=StepMagnitude

    elseif (time>=7.2) & (time<=7.4)
         fStep=StepMagnitude

    elseif (time>=8.4) & (time<=8.6)
         fStep=StepMagnitude

     elseif (time>=9.6) & (time<=9.8)
          fStep=StepMagnitude

    else
        fStep=0.0
        # k=0.0
    end
    return fStep
end


StepMagnitude=-600 #N or about 60 kg
tArb=0:0.001:10
fArb=RectangularStepLoading.(StepMagnitude,tArb)

plot(tArb,fArb/1700,title="step function",label=false,xaxis="time (s)",yaxis="p(t)") #plotting my step function

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)


#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m_person, m_trampoline, E, A, width, c, spl  = p

    #interpolate abitrary force
    ptArb = spl(t[1])

    if ptArb<=0.0
        m=m_person+m_trampoline
    else
        m=m_trampoline
    end

    #k=8E*A/width^3*u[1]^2
    k=2*m/(L*tan(θ))*u[1]^2

    ddu[1] = -c/m*du[1] - k/m*u[1] + ptArb/m

end


#solver controls
dt=0.01 #seconds
totalt=10 #seconds

#define initial conditions
uo=0.0  #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m_person, m_trampoline, E, A, width, c, spl))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)
t=sol.t

#plot results
plot!(t,u,linewidth=1,
    xaxis="t [sec.]",yaxis="u [meters]", label="u vs t",title="displacement graph")
