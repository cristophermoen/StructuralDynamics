
#' # Lecture 12 - Earthquake ground motion simulation

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

using DifferentialEquations

#' Earthquake ground motions are interesting because every earthquake around the world is different, both in magnitude and in time. 

#' Let's consider an earthquake ground motion from an earthquake in 1979 centered in Imperial Valley, California.  The ground accelerations are obtained from the [Center for Engineering Strong Motion Data](https://www.strongmotioncenter.org/cgi-bin/CESMD/iqr_dist_DM2.pl?IQRID=ImperialValley79&SFlag=0&Flag=3).  

#' Read in ground acceleration data, note dims=(rows,columns) in text file that has 711 lines and 8 columns.
using DelimitedFiles
earthquake_data = readdlm("/home/jrun/data/code/Imperial_Valley_1979.txt",Float64,dims=(711,8))

#' Rearrange ground accleration data as a vector.
earthquake_data=transpose(earthquake_data)   #transpose matrix
utt_g=reshape(earthquake_data,(5688,1))  #reshape matrix into a vector, note 8*711=5688
utt_g=vec(utt_g)   #tell Julia to make the data a 1D vector

#' Units come in as cm/sec^2 from strongmotion.org, let's change them to m/sec^2.
utt_g=utt_g ./ 100.0

#' Define time range of the earthquake assuming the first acceleration reading occurs at t=0.0 seconds.  The ground acceleration is provided every 0.010 seconds.
t_eq=collect(range(0,length=5688,stop=5687*0.01))   #total time is 5687*0.01 seconds

#' Plot ground acceleration over time.
using Plots
plot(t_eq,utt_g,linewidth=1,title="Imperial Valley 1979 earthquake",
    xaxis="time, t (sec.)",yaxis="ground acc., m/sec^2",legend=false)

#' Convert the earthquake ground motion from a discrete data set to a continous function.  
using Dierckx
earthquake = Spline1D(t_eq,utt_g)

#' Package up the important physical parameters.
p = [k, m, c, earthquake];

#' Now let's write the equation of motion.
function SpringMassDamperEOM(utt,ut,u,p,t)

    m, k, c, earthquake  = p

    utt_g = earthquake(t[1])      #ground acceleration at time t

    utt[1] = -c/m*ut[1] - k*u[1]/m -utt_g

end


#' Run an SDOF dynamic simulation with a specific natural period of vibration.
Tn=0.5  #seconds

#' Define a mass.
m=300    #kg

#' Calculate the stiffness to arrive at the Tn above.
fn=1/Tn
ωn=2*pi*fn   #radians/second
k=ωn^2*m    #N/m

#' Define a viscous damping ratio.
ζ=0.05      #viscous damping ratio
ccr=2*m*ωn  #critical viscous damping coefficient
c=ζ*ccr #viscous damping coefficient, units of kg/sec.

#' Define some solver controls.
dt=0.01 #seconds
totalt=60 #seconds

#' Define initial conditions.  They are not needed here.
u_o=[0.0] #meters
ut_o=[0.0] #m/sec.

#' Solve.
prob = SecondOrderODEProblem(SpringMassDamperEOM, ut_o, u_o, (0.,totalt), (m, k, c, earthquake))
solution = solve(prob, DPRKN8(),tstops=0:dt:totalt)

#' Plot.
u = (x->x[2]).(solution.u)  #displacement
ut = (x->x[1]).(solution.u)  #velocity
t = solution.t

plot(t,u, xlabel="time [sec.]", ylabel="displacement, u [m]", legend=false)


#' Create animation.

using CairoMakie

#' Define initial position for mass-spring animation.
u_position = Node((u[1],0))


#' Define window and plotting limits.
scene1=Scene(resolution=(1000,500))
Xmin=minimum(u)
Ymin=-1.0
Xmax=maximum(u)
Ymax=1.0
Xrange = Xmax - Xmin
Yrange = Ymax - Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#' Show the mass.
Makie.scatter!(scene1, lift(xy->Point2f[xy],u_position), show_axis = false, marker = [:circle], color= :red, markersize = 30.0)


#' Animate the mass.
record(scene1, "animation.mp4", range(1, stop = length(u), step=100)) do i

    if i == 1
        sleep(1)
    end

    u_position[]=(u[i],0)

    sleep(1/24)

end