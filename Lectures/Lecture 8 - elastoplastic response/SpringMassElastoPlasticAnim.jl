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
a=10000 #N
ω=0.2*ωn  #radians/sec.

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

StepMagnitude=0.0
StepTimeStart=0.0
StepTimeEnd=4.0
tArb=0:0.001:10
fArb=RectangularStepLoading.(StepMagnitude,StepTimeStart,StepTimeEnd,tArb)

#define abitrary force as a spline for later interpolation
spl = Spline1D(tArb,fArb)


#define elastoplastic response

uy=0.3  #yield displacement, meters

#set up global variables to track inelastic SDOF response
global up=0.0   #initialize plastic displacement, note global variable designation so that up is updated inside the SDOF function
global i=1   #counter for tracking fs
global fs_all=zeros(300000)   #define fs tracking vector, make vector really long to accommodate all the Julia time steps
global u_all=zeros(300000)   #define u tracking vector
global t_all=zeros(300000)   #define t tracking vector


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
          fs=k*(u[1]-up)    #calculate fs when displacement is elastic loading or unloading
    else
          fs=k*uy*sign(u[1]-up)   #calculate fs when displacement is plastic, on horizontal part of the curve
          global up=u[1]-uy*sign(u[1]-up)   #define how much plastic displacement occurs within a timestep and keep track of it globally
    end

      global fs_all[i]=fs
      global u_all[i]=u[1]
      global t_all[i]=t[1]
      global i=i+1

    ddu[1] = -c/m*du[1] - fs/m + ptForcing/m -F*SignF/m + ptArb/m
end

#solver controls
dt=0.01 #seconds
totalt=10 #seconds

#define initial conditions
uo=0.5*a/k  #meters
duo=ωn*a/k #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m, k, c, a, ω, F, spl,uy))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out u(t),du(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)   #angular velocity
t=sol.t

#clean zeros from intermediate solution
t_all=t_all[1:i-1]
u_all=u_all[1:i-1]
fs_all=fs_all[1:i-1]


#make some plots
using Makie


scene=Scene()
lines!(scene,t,u, color = :blue)
axis = scene[Axis]
axis[:names, :axisnames] = ("t [sec.]","u(t) [meters]")


scene=Scene()
lines!(t_all,fs_all, color = :blue)
axis = scene[Axis]
axis[:names, :axisnames] = ("t [sec.]","fs(t) [N]")

scene=Scene()
lines!(u_all,fs_all, color = :blue)
axis = scene[Axis]
axis[:names, :axisnames] = ("u [meters]","fs(t) [N]")


#animate u vs. fs plot with a red point follower

#tuple of x-y follower position, the Node command automatically updates later
#in the record function
XYPos=Node((u_all[1],fs_all[1]))

#define window and limits
scene=Scene(resolution=(1000,500))
Xmin=minimum(u_all)
Ymin=minimum(fs_all)
Xmax=maximum(u_all)
Ymax=maximum(fs_all)
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#plot the whole response
lines!(scene, u_all,fs_all, limits=limits, color= :blue)

#show pointer
scatter!(scene, lift(xy->Point2f0[xy],XYPos), marker = [:circle], limits=limits, color= :red, markersize=1000)

#animate
record(scene, "ufsAnimation.mp4", range(1, stop = length(u_all), step=1)) do i
     XYPos[]=(u_all[i],fs_all[i])
end
