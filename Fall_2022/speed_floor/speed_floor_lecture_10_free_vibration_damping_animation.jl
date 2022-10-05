using InstantFrame, DifferentialEquations


material = InstantFrame.Material(names=["steel"], E=[29500.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[0.0478*(36+24*2)], Iy=[37.1], Iz=[18.537], J=[0.001])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))


L = 16.0 * 12
num_elem = 10
x = range(0.0, L, num_elem+1)
numbers = range(1, num_elem+1)
coordinates = [(x[i], 0.0, 0.0) for i in eachindex(x)]
node = InstantFrame.Node(numbers=numbers, coordinates=coordinates)


element = InstantFrame.Element(numbers=1:num_elem, nodes=[(i, i+1) for i=1:num_elem], orientation=zeros(Float64, num_elem), connections=[("rigid", "rigid") for i=1:num_elem], cross_section=["beam" for i=1:num_elem], material=["steel" for i=1:num_elem])


#capital X is global X axis, capital Y is global Y axis, capital Z is global Z axis
#Inf means fixed, 0.0 means free 
nodes = 1:num_elem+1
uX = zeros(Float64, num_elem+1)
uX[1] = Inf

uY = zeros(Float64, num_elem+1)
uY[1] = Inf
uY[end] = Inf

uZ = [Inf for i=1:num_elem+1]

rX = [Inf for i=1:num_elem+1]  #now fixed for twist

rY = [Inf for i=1:num_elem+1]

rZ = zeros(Float64, num_elem+1)


support = InstantFrame.Support(nodes=nodes, stiffness=(uX=uX, uY=uY, uZ=uZ, rX=rX, rY=rY, rZ=rZ))

uniform_load = InstantFrame.UniformLoad(nothing)

point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "modal vibration")


#check natural frequency
#https://roymech.org/Useful_Tables/Vibrations/Natural_Vibrations.html
A = cross_section.A[1]
E = material.E[1]
I = cross_section.Iz[1]
L = 192.0
ρ = material.ρ[1]
m = ρ * A  #mass per unit length 

ωn = 9.87/L^2*sqrt(E*I/m)

#all good if true
isapprox(model.solution.ωn[1], ωn, rtol=0.05)



#define M, K, and C 

M = model.equations.M
K = model.equations.Ke

#' Calculate Raleigh damping C = ao*M + a1*K by assuming the modal damping ratio ζn for modes 1 and 2.

ζ1 = 0.05
ζ2 = 0.05

#' Find ao and a1.
a = (1/2 * [1/model.solution.ωn[1] model.solution.ωn[1]
            1/model.solution.ωn[2] model.solution.ωn[2]])  \ [ζ1; ζ2] 

# Define the Raleigh damping matrix.
# C = a[1] * M + a[2] * K   #a[2] is really small, which is causing problems
C = a[1] * M  #try just mass proportional damping

# Write the equation of motion.
function mdof(utt, ut, u, p, t)

    M, K, C = p

    utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:, 1]
    #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end


## Assume the initial displacement for free vibration is mode shape 1.

# Convert modal displaced shape into nodal displacements, easier to deal with.
mode_nodal_displacements = InstantFrame.define_nodal_displacements(node, model.solution.ϕ[1])

# Get global y displacement and check mode shape.
Y = [mode_nodal_displacements[i][2] for i in eachindex(mode_nodal_displacements)]
# plot(x, Y, markershape = :o)

# Normalize mode shape so that the maximum is unity.
mode_shape = model.solution.ϕ[1] ./ minimum(model.solution.ϕ[1])

#Check shape again.
mode_nodal_displacements = InstantFrame.define_nodal_displacements(node, mode_shape)
Y = [mode_nodal_displacements[i][2] for i in eachindex(mode_nodal_displacements)]
# plot(x, Y, markershape = :o)


#Define initial conditions.
u_0 = mode_shape
ut_0 = zeros(Float64, size(K,1)) 

# Define the time range for the simulation.
t_min = 0.0
t_max = 0.1
t_step = 0.01
t = range(t_min, t_max, step=t_step)

# Send in reduced matrices and vectors.
Mff = M[model.equations.free_dof, model.equations.free_dof]
Kff = K[model.equations.free_dof, model.equations.free_dof]
Cff = C[model.equations.free_dof, model.equations.free_dof]

p = (Mff, Kff, Cff)

u_0ff = u_0[model.equations.free_dof]
ut_0ff = ut_0[model.equations.free_dof]
                              
problem = SecondOrderODEProblem(mdof, ut_0ff, u_0ff, (t_min,t_max), p)

# Solve.
solution = solve(problem, DPRKN8(),tstops=t)


#Find all global-Y dof.

Y_dof = 2:6:length(node.numbers)*6

Y_dof_free = Y_dof[2:end-1]

node_index = [findfirst(num->num==Y_dof_free[i], model.equations.free_dof) for i in eachindex(Y_dof_free)]


# u_Y_t = (x->x[2]).(solution.u)

# u_Y = [(x->x[2][index]).(solution.u[i]) for i in eachindex(solution)]

# u_free = [(y->y[Y_dof_free[i]]).(solution.u) for i in eachindex(Y_dof_free)]

# #Plot solution
# plot(solution.t, u_midspan, legend=false)


#Work on animation

#Test out integrator, cool...
integ = init(problem, Tsit5())
# step!(integ, 0.01)
# Δ = (x->x[node_index]).(integ.u)

# Δ = get_beam_displacements(integ, node_index)


#Write up functions for animation

#get midspan displacement
# function get_midspan_displacement(integ, node_index)

#     Δ = integ.u[node_index]
#     return Δ

# end

#get all beam displacements
function get_beam_displacements(integ, node_index)

    Δ = integ.u.x[2][node_index]

    Δ = [0.0; Δ; 0.0]

    return Δ

end



#get displacement for a step in the solution
function progress_for_one_step!(integ, node_index)

    step!(integ, 0.01)
    Δ = get_beam_displacements(integ, node_index)

    return Δ

end


#initialize animation
using GLMakie
# node_index = findfirst(num->num==32, model.equations.free_dof)
Δ = u_0ff[node_index]

beam_points = [Observable(Point2f0(x[i], Δ[i])) for i in eachindex(Δ)]

fig = Figure(); display(fig)

ax = Axis(fig[1,1])

for i in eachindex(beam_points)
    GLMakie.scatter!(ax, beam_points[i]; marker=:circle, strokewidth=2, strokecolor=:red, color=:black, markersize = [8])
end

ax.title = "beam displacement"
ax.aspect = DataAspect()
GLMakie.ylims!(ax, -1.5, 1.5)


Δ = progress_for_one_step!(integ, node_index)


function animstep!(integ, beam_points, node_index)

    Δ = progress_for_one_step!(integ, node_index)
    beam_points[] = [Point2f0(x[i], Δ[i]) for i in eachindex(Δ)]

end

for i=1:10
    animstep!(integ, beam_points, node_index)
    sleep(0.001)
end


function make_animation(problem)

    integ = init(problem, Tsit5())
    node_index = findfirst(num->num==32, model.equations.free_dof)
    Δ = u_0ff[node_index]


    mid_point = Observable(Point2f0(0.0, Δ))

    fig = Figure(); display(fig)

    ax = Axis(fig[1,1])

    GLMakie.scatter!(ax, mid_point; marker=:circle, strokewidth=2, strokecolor=:red, color=:black, markersize = [8])

    ax.title = "midspan displacement"
    ax.aspect = DataAspect()
    GLMakie.ylims!(ax, -1.5, 1.5)

    return fig, integ, mid_point

end


fig, integ, mid_point = make_animation(problem)

frames = 1:200
record(fig, "speed_floor.mp4", frames; framerate = 60) do i
    for j = 1:5
        animstep!(integ, mid_point, node_index)
    end
end
