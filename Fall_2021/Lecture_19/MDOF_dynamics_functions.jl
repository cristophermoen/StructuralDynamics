using LinearAlgebra
using Dierckx
using DiffEqOperators
using DifferentialEquations

function calculate_member_length(member_definitions, node_geometry)

    L=zeros(length(member_definitions))

    for i=1:length(member_definitions)

        i_node = member_definitions[i][1]
        j_node = member_definitions[i][2]

        i_node_position = node_geometry[i_node,:]
        j_node_position = node_geometry[j_node,:]

        L[i] = norm(i_node_position-j_node_position)

    end

    return L

end


function define_member_property(member_definitions, member_property, property_order, property_type)

    A=zeros(length(member_definitions))

    for i = 1:length(member_definitions)

        A[i] = member_property[member_definitions[i][property_order]][property_type]

    end

    return A

end


function calculate_member_orientation(member_definitions, node_geometry)

    α=zeros(length(member_definitions))

    for i=1:length(member_definitions)

        i_node = member_definitions[i][1]
        j_node = member_definitions[i][2]

        i_node_position = node_geometry[i_node,:]
        j_node_position = node_geometry[j_node,:]

        member_vector = j_node_position - i_node_position
        horizontal_vector = [1; 0]

        α[i] =acos(dot(member_vector,horizontal_vector)/(norm(member_vector)*norm(horizontal_vector)))

    end

    return α

end

function define_local_Ke(A, E, Ix, L)

    k=zeros(6,6)

    k[1,1]=A*E/L
    k[1,4]=-A*E/L
    k[2,2]=12E*Ix/L^3
    k[2,3]=6E*Ix/L^2
    k[2,5]=-12E*Ix/L^3
    k[2,6]=6E*Ix/L^2
    k[3,3]=4E*Ix/L
    k[3,5]=-6E*Ix/L^2
    k[3,6]=2E*Ix/L
    k[4,4]=A*E/L
    k[5,5]=12*E*Ix/L^3
    k[5,6]=-6E*Ix/L^2
    k[6,6]=4E*Ix/L

    k[2:6,1]=k[1,2:6]
    k[3:6,2]=k[2,3:6]
    k[4:6,3]=k[3,4:6]
    k[5:6,4]=k[4,5:6]
    k[6,5]=k[5,6]

    return k

end


function define_rotation_matrix(α)

    θ=zeros(Float64, (6,6))

    θ[1,1]=cos(α)
    θ[1,2]=sin(α)
    θ[2,1]=-sin(α)
    θ[2,2]=cos(α)
    θ[3,3]=1.0

    θ[4:6,4:6]=θ[1:3,1:3]

    return θ

end

function transform_local_to_global_matrix(local_matrix, θ)

    global_matrix = transpose(θ) * local_matrix * θ

end


function define_local_Me(A,L,ρ)

    m=zeros(Float64, (6,6))

    m[1,1] = 140
    m[1,4] = 70
    m[2,2] = 156
    m[2,3] = 22L
    m[2,5] = 54
    m[2,6] = -13L
    m[3,3] = 4L^2
    m[3,5] = 13L
    m[3,6] = -3L^2
    m[4,4] = 140
    m[5,5] = 156
    m[5,6] = -22L
    m[6,6] = 4L^2

    m[2:6,1] = m[1,2:6]
    m[3:6,2] = m[2,3:6]
    m[4:6,3] = m[3,4:6]
    m[5:6,4] = m[4,5:6]
    m[6,5] = m[5,6]

    m =(ρ*A*L/420) * m

    return m

end

function define_fixed_dof(supports)

    global_fixed_dof=zeros(Int64,length(supports))

    for i=1:length(supports)

        node_number = supports[i][1]

        dof_start = 3 * (node_number - 1)

        global_fixed_dof[i] = dof_start + supports[i][2]

    end

    return global_fixed_dof

end


function assemble_global_matrix(node_geometry, member_definitions, global_element_matrix)

    num_nodes = size(node_geometry)[1]

    global_matrix = zeros(num_nodes * 3, num_nodes * 3)

    for i=1:length(member_definitions)

        i_node = member_definitions[i][1]
        j_node = member_definitions[i][2]

        i_node_dof = [1;2;3] .+(i_node - 1) * 3
        j_node_dof = [1;2;3] .+(j_node - 1) * 3

        global_element_dof = [i_node_dof; j_node_dof]

        global_matrix[global_element_dof, global_element_dof] = global_matrix[global_element_dof, global_element_dof] + global_element_matrix[i]

    end

    return global_matrix

end

function define_free_dof(supports)

    fixed_dof = define_fixed_dof(supports)

    all_dof=1:size(node_geometry)[1]*3

    free_dof=setdiff(all_dof,fixed_dof)

    return free_dof

end



function calculate_global_K(member_definitions, section_properties, material_properties, node_geometry, supports)

    #calculate member length
    L = calculate_member_length(member_definitions, node_geometry)

    #calculate member orientation
    α = calculate_member_orientation(member_definitions, node_geometry)

    #define section properties for stiffness matrix
    A = define_member_property(member_definitions, section_properties, 3, 1)
    E = define_member_property(member_definitions, material_properties, 4, 1)
    Ix = define_member_property(member_definitions, section_properties, 3, 2)

    #calculate local stiffness matrix for each element
    local_Ke = define_local_Ke.(A, E, Ix, L)

    #define rotation matrix for each element
    θ = define_rotation_matrix.(α)

    #rotate local stiffness matrix into global coordinates
    global_Ke = transform_local_to_global_matrix.(local_Ke, θ)

    #assemble global stiffness matrix
    global_K = assemble_global_matrix(node_geometry, member_definitions, global_Ke)

    #partition global stiffness matrix, only free DOF
    free_dof = define_free_dof(supports)
    K = global_K[free_dof, free_dof]

    return K

end


function calculate_global_M(member_definitions, section_properties, material_properties, node_geometry, supports)

    #define section and material properties
    A = define_member_property(member_definitions, section_properties, 3, 1)
    ρ = define_member_property(member_definitions, material_properties, 4, 2)

    #calculate member length
    L = calculate_member_length(member_definitions, node_geometry)

    #calculate member orientation
    α = calculate_member_orientation(member_definitions, node_geometry)

    #define local mass matrix
    local_Me = define_local_Me.(A,L,ρ)

    #define rotation matrix for each element
    θ = define_rotation_matrix.(α)

    #rotate local mass matrix into global coordinates
    global_Me = transform_local_to_global_matrix.(local_Me, θ)

    #assemble global mass matrix
    global_M = assemble_global_matrix(node_geometry, member_definitions, global_Me)

    #partition global mass matrix, only free DOF
    free_dof = define_free_dof(supports)
    M = global_M[free_dof, free_dof]

    return M

end