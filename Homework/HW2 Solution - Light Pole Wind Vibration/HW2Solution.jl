
#Assignment
#https://docs.google.com/document/d/1tgflJUNkUJc0_0swLOhlgYpnDIw_B6UfHjFKDTVEWwU/edit?usp=sharing

using OrdinaryDiffEq  #define Julia ODE package that we will use as solver
using Distributions   #for Gaussian wind distribution
using Dierckx  #this is for interpolating an arbitrary loading

#define constants
m=100    #kg
g=9.8  #m/sec^2
fn=51/36  #cycles per second
Tn=1/fn   #seconds
ωn=2*pi*fn   #radians/second


k=ωn^2*m    #N/m

ζ=0.01      #viscous damping ratio
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

#*****************
#random wind loading

#first quantify the mean wind gust frequency, the forcing frequency in this case
σGustTn=5.33 #seconds
σfGust=1/σGustTn #cycles per seconds

#mean and standard deviation of gust frequency
σωGust=σfGust*2*pi  #radians per second
μωGust=0.1*σωGust   #radians per second

#typically gust frequency is modeled as a Gaussian random variable
ωGustDist=Normal(σωGust,μωGust)

#now consider magnitude of the wind
#typical model  wind speed (t) = mean wind speed + gust(t)
#mean wind speed during a storm
vMean=60.0 #mph
#gust wind speed during a storm
vGustMean=30.0 #mph

#area of light pole
A=1*15  #ft^2

function WindLoadModel(vMean,vGustMean,ωGustDist,A,t)

    ωGust=rand(ωGustDist,1)

    #velocity as a function of time
    v=vMean+vGustMean*sin(ωGust[1]*t)

    #convert velocity to pressure
    pressure=0.00256*v^2    #psf

    #calculate force on SDOF system, assume constant pressure up the pole
    fWind=pressure*A  #lbs

    if fWind<0
        fWind=0
    end

    return fWind

end


tArb=0:0.01:50
fArb=WindLoadModel.(vMean,vGustMean,ωGustDist,A,tArb)

plot(tArb,fArb)

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)


#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, a, ω, F, spl  = p

    #harmonic forcing function
    ptForcing=a*sin(ω*t[1])

    #interpolate abitrary force
    ptArb = spl(t[1])    #interpolate abitrary force

    #friction damping
    if du[1]>0
        SignF=1
    else
        SignF=-1
    end

    ddu[1] = -c/m*du[1] - k/m*u[1] + ptForcing/m -F*SignF/m + ptArb/m
end

#solver controls
dt=0.001 #seconds
totalt=50 #seconds

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
