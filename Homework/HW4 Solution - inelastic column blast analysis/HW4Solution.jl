#define Julia ODE package that we will use as solver
using OrdinaryDiffEq
using Dierckx  #this is for interpolating an aribtrary loading

#define mass and stiffness
m=460    #kg, 15 ft. tall W14x61 and 0.1m facade
g=9.8  #m/sec^2
k=42800988    #N/m for SDOF approximation of column, k=384EI/(5L^3)

#calculate frequency and period
ωn=sqrt(k/m)  #radians per second
fn=ωn/(2*pi)  #cycles per second
Tn=1/fn       #period, seconds

ζ=0.05      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

μ=0.00  #friction coefficient
N=m*g   #normal force
F=μ*N   #friction force

#forcing function
a=0 #N
ω=0.2*ωn  #radians/sec.

#define arbitrary function
#in this case a blast load
using DelimitedFiles

# assume txt file to import is in the same folder as this Julia script
cd(@__DIR__)  #change working directory
fArb= vec(readdlm("blastforce.txt"))  #load the blast force, psi
tArb=collect(range(0,length=5001,stop=0.5))  #define time range

TribArea=15*40*12^2  #in^2, facade tributary area

fArb=fArb*TribArea*4.448   #convert psi to lbs, then lbs to N

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)

#define elastoplastic response
uy=0.02127  #yield displacement, meters, calculated just from bending

#set up global variables to track inelastic SDOF response
global up=0.0   #initialize plastic displacement, note global variable designation so that up is updated inside the SDOF function
global i=1   #counter for tracking fs
global fS_all=zeros(300000)   #define fs tracking vector, make vector really long to accommodate all the Julia time steps
global u_all=zeros(300000)   #define u tracking vector
global du_all=zeros(300000)   #define velocity tracking vector
global ddu_all=zeros(300000)   #define acceleration tracking vector
global t_all=zeros(300000)   #define t tracking vector
global up_all=zeros(300000)   #define up tracking vector


#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, a, ω, F, spl, uy  = p

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

    #define elastoplastic response
    if abs(u[1]-up) < uy
          fS=k*(u[1]-up)    #calculate fs when displacement is elastic loading or unloading
    else
          fS=k*uy*sign(u[1]-up)   #calculate fs when displacement is plastic, on horizontal part of the curve
          global up=u[1]-uy*sign(u[1]-up)   #define how much plastic displacement occurs within a timestep and keep track of it globally
          #up can be positive or negative and it is calculated cumulatively
          #plastic deformation to the right is positive, deformation to the left is negative

    end

      global fS_all[i]=fS
      global u_all[i]=u[1]
      global du_all[i]=du[1]
      global t_all[i]=t[1]
      global up_all[i]=up


    ddu[1] = -c/m*du[1] - fS/m + ptForcing/m -F*SignF/m + ptArb/m

    global ddu_all[i]=ddu[1]
    global i=i+1

end

#solver controls
dt=0.001 #seconds
totalt=1 #seconds

#define initial conditions
uo=0.0  #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m, k, c, a, ω, F, spl,uy))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)   #angular velocity
t=sol.t

#clean zeros from intermediate solution
t=t_all[1:i-1]
u=u_all[1:i-1]
du=du_all[1:i-1]
ddu=ddu_all[1:i-1]
fS=fS_all[1:i-1]
up_all=up_all[1:i-1]

#plot displacement vs. time
#note plastic deformation of about 0.5 meters
using Plots
plot(t,u,xaxis="t [sec.]",yaxis="u [m]", legend=false)

#calculate steps in u along time history
Δu=[0; diff(u)]

#calculate inertial force
fI=m*ddu
#calculate viscous damping force
fD=c*du

#and fS has already been calculated, remember it changes over the time history


#plot D'Alembert's forces vs. time
p1=plot(t,fI,xaxis="t (sec.)",yaxis="fI, N")
p2=plot(t,fS,xaxis="t (sec.)",yaxis="fS, N")
p3=plot(t,fD, xaxis="t (sec.)",yaxis="fD, N")
plot(p1,p2,p3,layout=(3,1),legend=false)


#calculate energy quantities

#write a little function to help with this
#integrate little chunks of P vs. Δu
function CalculateEnergy(P,Δu)
    energy=zeros(length(u))
    for i=2:length(P)
        energy[i]=P[i]*Δu[i]
    end
    return energy
end

#chunks of energy
dEnergyI=CalculateEnergy(fI,Δu)  #inertial
dEnergyS=CalculateEnergy(fS,Δu)  #spring
dEnergyD=CalculateEnergy(fD,Δu)  #viscous damping

#cumulative energy
CumulativeEnergyI=accumulate(+, dEnergyI)
CumulativeEnergyS=accumulate(+, dEnergyS)
CumulativeEnergyD=accumulate(+, dEnergyD)


#compare inertial energy to something we know
#E=1/2*m*du(0)^2

#calculate initial velocity from blast
plot(t,du)
duMax=maximum(du)  #m/s^2,

#compare to inertial energy calculated from m*ddu
plot(t,CumulativeEnergyI,xaxis="t (sec.)",yaxis="CeI, N-m")
KineticEnergy=1/2*m*duMax^2  #N-m
plot!([0; 1],[KineticEnergy; KineticEnergy],linecolor = :red)


#plot cumulative energy vs. time
p1=plot(t,CumulativeEnergyI,xaxis="t (sec.)",yaxis="CeI, N-m")
p2=plot(t,CumulativeEnergyS,xaxis="t (sec.)",yaxis="CeS, N-m")
p3=plot(t,CumulativeEnergyD,xaxis="t (sec.)",yaxis="CeD, N-m")
plot(p1,p2,p3,layout=(3,1),legend=false)
