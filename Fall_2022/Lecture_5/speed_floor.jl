using InstantFrame, DifferentialEquations


material = InstantFrame.Material(names=["steel"], E=[29500.0], ν=[0.3], ρ=[492.0 / 32.17 / 12^4 / 1000.0])  ##ρ = kilo-lbs * s^2 / in^4

cross_section = InstantFrame.CrossSection(names=["beam"], A=[9.12], Iy=[37.1], Iz=[110.0], J=[0.001])

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

rX = zeros(Float64, num_elem+1)
rX[1] = Inf

rY = [Inf for i=1:num_elem+1]

rZ = zeros(Float64, num_elem+1)


support = InstantFrame.Support(nodes=nodes, stiffness=(uX=uX, uY=uY, uZ=uZ, rX=rX, rY=rY, rZ=rZ))



uniform_load = InstantFrame.UniformLoad(nothing)
point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "modal vibration")

#Cassie - here are [K] and [M]
model.equations.Ke
model.equations.M 

#The global free dof are also provide based on the boundary conditions I assigned above.
model.equations.free_dof 

#which means we can get the partitioned matrix to put into the ODE solver
Ke_ff = model.equations.Ke[model.equations.free_dof, model.equations.free_dof]
M_ff = model.equations.M[model.equations.free_dof, model.equations.free_dof]

#The eigenvalue problem is solved as well, which means we get ωn and ϕn for each mode.   

model.solution.ωn
model.solution.ϕ

# We should check ωn,1 against an analytical solution, just to make sure InstantFrame is working. 

# All for now.   I will work on animation for next Tuesday since we will need that.
