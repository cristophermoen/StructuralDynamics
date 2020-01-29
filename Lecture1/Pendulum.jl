
#define Julia ODE package that we will use as solver
using OrdinaryDiffEq

#define constants
g=9.8  #m/s^2
L=0.5 #meters

#write equation of motion for a pendulum
function PendulumEOM(ddθ,dθ,θ,p,t)
    g, L = p
    ddθ[1] = -g*sin(θ[1])/L
end

#solver controls
dt=0.01 #seconds
totalt=10 #seconds

#define initial conditions
θo=pi/12
dθo=0.0

#solve the ODE
prob = SecondOrderODEProblem(PendulumEOM, [dθo], [θo], (0.,totalt), (g,L))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out θ(t),dθ(t), and t from solution
dθ=first.(sol.u)
θ=last.(sol.u)   #angular velocity
t=sol.t

#plot results
using Plots    #start Julia plot package
plot(t,θ,linewidth=1,title="pendulum free vibration",
    xaxis="time, t (sec.)",yaxis="rotation, \\theta (t) (radians)", legend=false)
