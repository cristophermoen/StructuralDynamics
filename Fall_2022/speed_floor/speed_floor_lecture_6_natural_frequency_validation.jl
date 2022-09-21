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

model.solution.ωn[1]/(2π)

#define M, K, and C 

M = model.equations.M
K = model.equations.Ke

#' Calculate Raleigh damping C = ao*M + a1*K by assuming the modal damping ratio ζn for modes 1 and 3.

ζ1 = 0.05
ζ3 = 0.05

#' Find ao and a1.
a = (1/2 * [1/model.solution.ωn[1] model.solution.ωn[1]
            1/model.solution.ωn[3] model.solution.ωn[3]])  \ [ζ1; ζ3] 

# Define the Raleigh damping matrix.
C = a[1] * M + a[2] * K


# Write the equation of motion.
function mdof(utt, ut, u, p, t)

    # M, K, C = p
    M, K = p

    # utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:, 1]
    utt[:,1] = - M \ K * u[:,1]
    #https://stackoverflow.com/questions/54890422/inv-versus-on-julia  See this link for discussion of inv(M) * K  vs. M \ K.  

end


# Assume the initial displacement for free vibration is an initial velocity at midspan.
u_0 = zeros(Float64, size(K,1))  
ut_0 = zeros(Float64, size(K,1)) 

midspan_Y_dof = 5*6 + 2
# ut_0[midspan_Y_dof] = -1.0
u_0[midspan_Y_dof] = -1.0

# Define the time range for the simulation.
t = 0.0:0.01:0.1

# Send in reduced matrices and vectors.
Mff = M[model.equations.free_dof, model.equations.free_dof]
Kff = K[model.equations.free_dof, model.equations.free_dof]
Cff = C[model.equations.free_dof, model.equations.free_dof]

p = (Mff, Kff)

u_0ff = u_0[model.equations.free_dof]
ut_0ff = ut_0[model.equations.free_dof]
                              
problem = SecondOrderODEProblem(mdof, u_0ff, ut_0ff, (0.,0.1), p)

#' Solve.
solution = solve(problem, DPRKN8(),tstops=t)

index = findfirst(num->num==32, model.equations.free_dof)

length(model.equations.free_dof)

u_midspan=(x->x[19]).(solution.u)
t_plot=solution.t

using Plots
plot(t_plot, u_midspan)