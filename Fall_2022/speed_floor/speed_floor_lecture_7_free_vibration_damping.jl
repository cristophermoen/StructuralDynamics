using InstantFrame, DifferentialEquations, Plots


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
plot(x, Y, markershape = :o)

# Normalize mode shape so that the maximum is unity.
mode_shape = model.solution.ϕ[1] ./ minimum(model.solution.ϕ[1])

#Check shape again.
mode_nodal_displacements = InstantFrame.define_nodal_displacements(node, mode_shape)
Y = [mode_nodal_displacements[i][2] for i in eachindex(mode_nodal_displacements)]
plot(x, Y, markershape = :o)


#Define initial conditions.
u_0 = mode_shape
ut_0 = zeros(Float64, size(K,1)) 

# Define the time range for the simulation.
t_min = 0.0
t_max = 0.5
t_step = 0.1
t = range(t_min, t_max, step=t_step)

# Send in reduced matrices and vectors.
Mff = M[model.equations.free_dof, model.equations.free_dof]
Kff = K[model.equations.free_dof, model.equations.free_dof]
Cff = C[model.equations.free_dof, model.equations.free_dof]

p = (Mff, Kff, Cff)

u_0ff = u_0[model.equations.free_dof]
ut_0ff = ut_0[model.equations.free_dof]
                              
problem = SecondOrderODEProblem(mdof, u_0ff, ut_0ff, (t_min,t_max), p)

# Solve.
solution = solve(problem, DPRKN8(),tstops=t)


#Find midspan global-Y dof.
index = findfirst(num->num==32, model.equations.free_dof)
u_midspan=(x->x[index]).(solution.u)

#Plot solution
plot(solution.t, u_midspan, legend=false)


