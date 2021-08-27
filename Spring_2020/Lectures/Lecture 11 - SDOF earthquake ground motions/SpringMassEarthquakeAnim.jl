using DelimitedFiles #use for importing earthquake acceleration file
using OrdinaryDiffEq #define Julia ODE package that we will use as solver
using Dierckx        #this is for interpolating an aribtrary loading


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


#make some plots
using Makie

#define initial position for mass-spring animation
XYPos1=Node((u[1],0))

#define initial position for displacement vs. time animation
XYPos2=Node((t[1],u[1]))

#define window and limits
scene1=Scene(resolution=(1000,500))
Xmin=minimum(u)
Ymin=-1
Xmax=maximum(u)
Ymax=1
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#show mass
scatter!(scene1, lift(xy->Point2f0[xy],XYPos1), show_axis = false, marker = [:rect], limits=limits, color= :red, markersize=0.5)


#define window and limits
scene2=Scene(resolution=(1000,500))
Xmin=minimum(t)
Ymin=minimum(u)
Xmax=maximum(t)
Ymax=maximum(u)
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#plot the whole response
lines!(scene2, t,u, limits=limits, color= :blue)

axis=scene2[Axis]

axis[:names, :axisnames] = ("time [sec.]", "u [m]")

#show pointer
scatter!(scene2, lift(xy->Point2f0[xy],XYPos2), marker = [:circle], limits=limits, color= :red, markersize=2)


#group the scenes
scene3=vbox(scene2,scene1)

#animate
#adjust step in loop to match actual time with animation time 
record(scene3, "animation.mp4", range(1, stop = length(u), step=4)) do i

    val, looptime, bytes, gctime, memallocs=@timed begin  #measure time for each loop
    XYPos1[]=(u[i],0)
    XYPos2[]=(t[i],u[i])
    end

    # sleep(dt-looptime)  #hold to correct loop time to 'real' time in animation

end
