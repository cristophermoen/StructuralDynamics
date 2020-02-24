
#Homework 3 - column blast loading
#https://docs.google.com/document/d/1DRJJr9Iu4DmrNjPNUW9sqp_prFWkVgbdZNP_BMTM_gc/edit?usp=sharing


#define Julia ODE package that we will use as solver
using OrdinaryDiffEq
using Dierckx  #this is for interpolating an aribtrary loading

#define constants
m=460    #kg
g=9.8  #m/sec^2

#assume blast force concentrated at middle of simply supported column
#k=48EI/L^3   concentrated load a center of a beam
k=27000000  #kN/m

ωn=sqrt(k/m)
fn=wn/(2*pi)
Tn=1/fn

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
#in this case a blast load
using DelimitedFiles

# assume txt file to import is in the same folder as this Julia script
cd(@__DIR__)  #change working directory
fArb= vec(readdlm("blastforce.txt"))  #load the blast force, psi
tArb=collect(range(0,length=5001,stop=0.5))  #define time range

plot(tArb,fArb,xlims=(0,0.05))

TribArea=15*40*12^2  #in^2, facade tributary area

td=0.007 #seconds, approximate duration of blast loading
tdNorm=td/Tn  #compare this to Chopra Figure 4.7.2, it is short

fArb=fArb*TribArea*4.448   #convert psi to lbs, then lbs to N

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)


#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, a, ω, F, spl  = p

    #harmonic forcing function
    ptForcing=a*sin(ω*t[1])

    #interpolate abitrary force, send it to ODE solver
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

#calculate and plot static displacement
po=maximum(fArb)
usto=po/k
plot!([0; 10],[usto; usto], linecolor= :red)
