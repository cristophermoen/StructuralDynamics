using DelimitedFiles #use for importing earthquake acceleration file
using OrdinaryDiffEq #define Julia ODE package that we will use as solver
using Dierckx        #this is for interpolating an aribtrary loading
using Plots

cd(@__DIR__)  #change working directory to that of this file

#read in ground acceleration data, note dims=(rows,columns) in text file that has 711 lines and 8 columns
#Imperial Valley 1979, Imperial County Center Grounds ground motion
data = readdlm("IV1979acc.txt",Float64,dims=(711,8))

#rearrange ground accleration data as a vector
data=transpose(data)   #transpose matrix
ddu_g=reshape(data,(5688,1))  #reshape matrix into a vector, note 8*711=5688
ddu_g=vec(ddu_g)   #tell Julia to make the data a 1D vector

#units come in as cm/sec^2 from strongmotion.org, let's change them to m/sec^2
ddu_g=ddu_g/100

#define time range of the earthquake assuming the first acc reading occurs at t=0 seconds
t_eq=collect(range(0,length=5688,stop=5687*0.01))   #total time is 5687*0.01 seconds

#plot ground acceleration over time
plot(t_eq,ddu_g,linewidth=1,title="Imperial Valley 1979 earthquake",
    xaxis="time, t (sec.)",yaxis="ground acc., m/sec^2",legend=false)

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

#calculate frequency and period
ωn=sqrt(k/m)  #radians per second
fn=ωn/(2*pi)  #cycles per second
Tn=1/fn       #period, seconds

#define ground acceleration as a spline for later interpolation
spl = Spline1D(t_eq,ddu_g)

#write equation of motion
function SpringMassDamperEOM(ddu,du,u,p,t)
    m, k, c, spl  = p

    ddu_g = spl(t[1])      #ground acceleration at time t

    ddu[1] = -c/m*du[1] - k*u[1]/m -ddu_g

end

#solver controls
dt=0.01 #seconds
totalt=60 #seconds

#define initial conditions
uo=0.0 #meters
duo=0.0 #m/sec.

#solve the ODE
prob = SecondOrderODEProblem(SpringMassDamperEOM, [duo], [uo], (0.,totalt), (m, k, c, spl))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out du(t), u(t), and t from solution
du=first.(sol.u)
u=last.(sol.u)
t=sol.t

#plot displacement vs. time
plot(t,u, xlabel="time [sec.]", ylabel="displacement, u [m]", legend=false)
