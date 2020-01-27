using DifferentialEquations

g=9.8  #m/s^2
L=0.1 #meters


#write equation of motion for a pendulum
function PendulumEOM(ddθ,θ,p,t)
    g = p
    ddθ[1] = -g[1]*sin(θ[1])/L
end
 #                                           θ0       dθ0   trange
  prob = SecondOrderODEProblem(PendulumEOM, [pi/12.], [0.], (0.,2.),g)
  sol = solve(prob, DPRKN6())

#separate out u(t),u_dot(t), and t from sol tuple
u=first.(sol.u)
u_dot=last.(sol.u)
t=sol.t

using Plots    #start Julia plot package

#plot displacement over time
plot(t,u,linewidth=1,title="SDOF solution with random excitation",
    xaxis="time, t (sec.)",yaxis="diplacement, u(t) (m)",label="My Thin Line!") # legend=false
