
#define Julia ODE package that we will use as solver
using OrdinaryDiffEq

#define constants
m=0.035    #kg
g=9.8  #m/s^2
L=0.5 #meters
c=0.018  #wind resistance parameter

#write equation of motion for a pendulum
function PendulumEOM(ddθ,dθ,θ,p,t)
    g, L, c, m = p
    ddθ[1] = -g*sin(θ[1])/L - c*dθ[1]/m
end

#solver controls
dt=0.01 #seconds
totalt=20 #seconds

#define initial conditions
θo=pi/4
dθo=0.0

#solve the ODE
prob = SecondOrderODEProblem(PendulumEOM, [dθo], [θo], (0.,totalt), (g,L,c, m))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out θ(t),dθ(t), and t from solution
dθ=first.(sol.u) #angular velocity
θ=last.(sol.u)
t=sol.t




#plot results
# using Plots    #start Julia plot package
# plot(t,θ,linewidth=1,title="pendulum free vibration",
#     xaxis="time, t (sec.)",yaxis="rotation, \\theta (t) (radians)", legend=false)




using Makie

 scene = Scene()
 mytheta = Node(0.0)
 xy(θt,L)=(L*sin(θt),L*cos(θt))
 scene = lines!(
     scene,
     lift(t -> xy.(range(0, stop = 2pi, length = 50), t), mytime),
     color = :blue)
 p1 = scene[end];
 N = 100
 record(scene, "output.mp4", range(0, stop = 4pi, length = N)) do i
     mytime[] = i
end





# using Makie
# scene = Scene()
#
# xy(θt,L)=(L*sin(θt),L*cos(θt))
#
# theta_node=Node(θ[1])
#
# cdm=lift(θt-> xy(θt, L), theta_node)
# points=[(0.0,0.0),cdm.val]
#
# p1=scatter!(scene,points)
#
# push!(theta_node,θ[1000])
#
# p1=scatter!(scene,[0.0],[0.0])[end]
# p2=scatter!(scene, lift(θt-> xy(θt, L), theta_node))[end]
#
# record(scene, "output.mp4", range(1, stop = 10, length = 10)) do i
#      push!(theta_node, θ[i])
# end

#
# Updateθ = Node(θ[1])
#
# push!(Updateθ,θ[200])
# y=lift(θt-> xy.(θt, L), Updateθ)



# xy1=(0.0,0.0)
# xy2=lift(θt-> xy.(θt, L), Updateθ)
#
# cdm=Point2f0(lift(θt-> xy.(θt, L), Updateθ))
#
# points = [
#     Point2f0(xy1[1], xy1[2]) => Point2f0(xy2[1], xy2[2]);
#     ]
#
# p1=scatter!(scene, [xy1[1]], [xy1[2]])[end]
# cdm=lift(θt-> xy.(θt, L), sendθ)
# # p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]])[end]
# p2=scatter!(scene, lift(θt-> xy.(θt, L), sendθ))
# )[end]
# points = lift(p1[1], p2[1]) do pos1, pos2
#      map((a, b)-> (a, b), pos1, pos2)
# end
# linesegments!(scene, points)
#
# i=1000
# push!(time_node, t[i])
# p1=scatter!(scene, [xy1[1]], [xy1[2]])[end]
# p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]])[end]
#
#
#
# points = lift(p1[1], p2[1]) do pos1, pos2
#      map((a, b)-> (a, b), pos1, pos2)
# end
# linesegments!(scene, points)
#
# scene = Scene()
# record(scene, "output.mp4", 1:length(t)) do i
#      push!(time_node, t[i])
#      # p1=scatter!(scene, [xy1[1]], [xy1[2]],show_axis=false)[end]
#      # p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]],show_axis=false)[end]
#      # points = lift(p1[1], p2[1]) do pos1, pos2
#      #      map((a, b)-> (a, b), pos1, pos2)
#      # end
#      # linesegments!(scene, points)
#      # scene=Scene()
# end



# xy(θt,L)=(L*sin(θt),L*cos(θt))
# i=1
# time_node = Node(t[i])
# xy1=(0.0,0.0)
# xy2=lift(θt-> xy.(θt, L), time_node)
#
#
#
# p1=scatter!(scene, [xy1[1]], [xy1[2]])[end]
# p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]])[end]
# points = lift(p1[1], p2[1]) do pos1, pos2
#      map((a, b)-> (a, b), pos1, pos2)
# end
# linesegments!(scene, points)
#
# i=1000
# push!(time_node, t[i])
# p1=scatter!(scene, [xy1[1]], [xy1[2]])[end]
# p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]])[end]
# points = lift(p1[1], p2[1]) do pos1, pos2
#      map((a, b)-> (a, b), pos1, pos2)
# end
# linesegments!(scene, points)
#
# scene = Scene()
# record(scene, "output.mp4", 1:length(t)) do i
#      push!(time_node, t[i])
#      # p1=scatter!(scene, [xy1[1]], [xy1[2]],show_axis=false)[end]
#      # p2=scatter!(scene, [xy2.val[1]],[xy2.val[2]],show_axis=false)[end]
#      # points = lift(p1[1], p2[1]) do pos1, pos2
#      #      map((a, b)-> (a, b), pos1, pos2)
#      # end
#      # linesegments!(scene, points)
#      # scene=Scene()
# end

# linesegments!(scene, points)
# y1=lift(θt-> y.(θt, L), time_node)

# cdm=lift(t-> f.(t, range(0, stop = 2pi, length = 50), 1), time_node)


# scene = Scene()
#  f(t, v, s) = (sin(v + t) * s, cos(v + t) * s)
#  time_node = Node(0.0)
#  p1 = scatter!(scene, lift(t-> f.(t, range(0, stop = 2pi, length = 50), 1), time_node))[end]
#  p2 = scatter!(scene, lift(t-> f.(t * 2.0, range(0, stop = 2pi, length = 50), 1.5), time_node))[end]
#  points = lift(p1[1], p2[1]) do pos1, pos2
#      map((a, b)-> (a, b), pos1, pos2)
# end
#  linesegments!(scene, points)
#  N = 150
#  record(scene, "output.mp4", range(0, stop = 10, length = N)) do i
#      push!(time_node, i)
# end
