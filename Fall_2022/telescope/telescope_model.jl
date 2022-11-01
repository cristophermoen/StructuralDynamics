using MAT, InstantFrame, GLMakie

telescope_mat_file = "/Users/crismoen/Documents/dev/StructuralDynamics/Fall_2022/telescope/T6.mat"

telescope = matread(telescope_mat_file)

show(keys(telescope))


telescope["node_info"]
telescope["elem_info"]
telescope["sect_info"]

telescope["fixity_info"]


#7900 kg/m^3
#7.9E-6 kg/mm^3

material = InstantFrame.Material(names=["steel"], E=[200.0], ν=[0.3], ρ=[7.9E-6])  #ρ = kg/mm^3

cross_section = InstantFrame.CrossSection(names=["hampPipe", "RHSS7625x328"], A=[6580.0, 4522.572], Iy=[45785456.0, 19604500.1646], Iz=[45785456.0, 19604500.1646], J=[91154682.0, 39167377.1866])

connection = InstantFrame.Connection(names=["rigid"], stiffness=(ux=[Inf], uy=[Inf], uz=[Inf], rx=[Inf], ry=[Inf], rz=[Inf]))


coordinates = [(telescope["node_info"][i,1], telescope["node_info"][i,3], telescope["node_info"][i,2]) for i=1:size(telescope["node_info"],1)]
node = InstantFrame.Node(numbers=1:length(coordinates), coordinates=coordinates)

num_elem = size(telescope["elem_info"],1)

element_cross_sections = Vector{String}(undef, num_elem)
for i=1:num_elem

    if telescope["elem_info"][i,3] == 1.0
        element_cross_sections[i] = "hampPipe"
    else
        element_cross_sections[i] = "RHSS7625x328"
    end

end

element_connectivity = [(Int(telescope["elem_info"][i,1]), Int(telescope["elem_info"][i,2])) for i=1:num_elem]
element = InstantFrame.Element(numbers=1:length(element_connectivity), nodes=element_connectivity, orientation=zeros(Float64, length(element_connectivity)), connections=[("rigid", "rigid") for i in eachindex(element_connectivity)], cross_section=element_cross_sections, material=["steel" for i in eachindex(element_connectivity)])

support = InstantFrame.Support(nodes=[1, 2, 3, 4], stiffness=(uX=[Inf,Inf,Inf,Inf], uY=[Inf,Inf, Inf,Inf], uZ=[Inf,Inf,Inf,Inf], rX=[Inf,Inf,Inf,Inf], rY=[Inf,Inf,Inf,Inf], rZ=[Inf,Inf,Inf,Inf]))

# uniform_load = InstantFrame.UniformLoad(labels=["test"], elements=[2], loads=(qX=zeros(Float64, 1), qY=zeros(Float64, 1), qZ=ones(Float64, 1)*-1000.0, mX=zeros(Float64, 1), mY=zeros(Float64, 1), mZ=zeros(Float64, 1)))
uniform_load = InstantFrame.UniformLoad(nothing)

# point_load = InstantFrame.PointLoad(labels = ["test"], nodes=[2], loads=(FX=[0.0], FY=[0.0], FZ=[0.0], MX=[0.0], MY=[0.0], MZ=[0.0]))
point_load = InstantFrame.PointLoad(nothing)

model = InstantFrame.solve(node, cross_section, material, connection, element, support, uniform_load, point_load, analysis_type = "first order")


###Visualize####

element_nodal_coords = InstantFrame.UI.define_element_nodal_start_end_coordinates(element, node)

X, Y, Z = InstantFrame.UI.get_node_XYZ(node)

X_range = abs(maximum(X) - minimum(X))
Y_range = abs(maximum(Y) - minimum(Y))
Z_range = abs(maximum(Z) - minimum(Z))



figure = Figure()
ax = Axis3(figure[1,1])
ax.aspect = (1.0, Y_range/X_range, Z_range/X_range)
ax.yticks = WilkinsonTicks(2)
ax.xticks = WilkinsonTicks(2)
# ylims!(ax, -50.0, 50.0)
# zlims!(ax, -50.0, 50.0)
# ylims!(ax, -50.0, 50.0)
# zlims!(ax, -50.0, 50.0)

color = :gray
InstantFrame.UI.show_elements!(ax, element_nodal_coords, color)
figure

markersize = 10
color = :blue
InstantFrame.UI.show_nodes!(ax, X, Y, Z, markersize, color)
figure

textsize = 10
color = :pink
InstantFrame.UI.show_element_numbers!(ax, element, node, textsize, color)
figure

textsize = 10
color = :black
InstantFrame.UI.show_node_numbers!(ax, node, textsize, color)
figure

# function show_element_numbers!(ax, element, node, textsize, color)

#     for i in eachindex(element.numbers)

#         start_index = findfirst(num->num==element.nodes[i][1], node.numbers)
#         end_index = findfirst(num->num==element.nodes[i][2], node.numbers)

#         Δ = node.coordinates[end_index] .- node.coordinates[start_index]

#         text_location = node.coordinates[start_index] .+ Δ./2

#         text!(ax,
#             [Point3f(text_location[1], text_location[2], text_location[3])],
#             text = [string(element.numbers[i])],
#             # rotation = [i / 7 * 1.5pi for i in 1:7],
#             color = color,
#             # align = (:left, :baseline),
#             textsize = textsize,
#             # markerspace = :data
#         )
#     end

# end











# function show_node_numbers!(ax, node, textsize, color)

#     text!(ax,
#         [Point3f(node.coordinates[i][1], node.coordinates[i][2], node.coordinates[i][3]) for i in eachindex(node.coordinates)],
#         text = [string(node.numbers[i]) for i in eachindex(node.numbers)],
#         # rotation = [i / 7 * 1.5pi for i in 1:7],
#         color = color,
#         # align = (:left, :baseline),
#         textsize = textsize,
#         # markerspace = :data
#     )

# end



ax.azimuth = 0.5π
ax.elevation = -0.5π
figure