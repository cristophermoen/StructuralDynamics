#' # Lecture 19 - Beam vibration using a consistent mass matrix

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' Simulate the vibration of a simply supported beam with an initial imposed
#displacement at midspan.   Use a 'stiffness method' approach where the beam is discretized into elements, the element stiffness and mass matrices are defined in local coordinates, then transformed to global coordinates, and assembled into global K and M matrices.   

#' Bring in functions to calculate global K and M.
include("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_19/MDOF_dynamics_functions.jl")

#' Define beam properties.

                      #A Ix b(depth)
section_properties = [(1.0, 1.0, 2.0),   #Section property 1
                     (2.0, 2.4, 220)]   #Section property 2

#E  ρ (unit weight)    N/m^2  #kg/m^3 in this problem
material_properties = [(1.0, 1.0),    #steel
                      (1E+5, 3000)]     #concrete

#x y
node_geometry = [0.0  0.0;   #Node 1
                12.0 0.0;   #Node 2
                24.0 0.0]   #Node 3

#iNode  jNode SectionProperties MaterialProperties
member_definitions = [(1, 2, 1, 1),   #Member 1
                     (2, 3, 1, 1)]   #Member 2

          # node number, DOF u=1 v=2 ϕ=3
supports = [(1,1),
            (1,2),
            (3,2)]

#' Write the equation of motion.
function mdof(utt, ut, u, p, t)

  M, K, C = p

  utt[:,1] = - M \ K * u[:,1] - M \ C * ut[:, 1]

end


K = calculate_global_K(member_definitions, section_properties, material_properties, node_geometry, supports)

M = calculate_global_M(member_definitions, section_properties, material_properties, node_geometry, supports)

#' Define the damping matrix $C$ using Raleigh damping.

#' Calculate natural frequencies.
ωn_squared=eigvals(K, M)
ωn=sqrt.(ωn_squared)

#' Define modal damping ratios for Raleigh damping.
ωi=ωn[1]   #set first frequency anchor point i
ωj=ωn[3]   #set second frequency anchor point j
ζi=0.05    #modal viscous damping ratio in mode i
ζj=0.05    #modal viscous damping ratio in mode j

a0,a1=2*inv([1/ωi ωi;1/ωj ωj])*[ζi;ζj]  #Eq. 11.4.9 from Chopra
C=a0*M+a1*K

 #                                    u_dot0                       u_0                         trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.; 0.; 0.; 0.], [0.0; 0.0; 0.1; 0.0; 0.0; 0.0], (0.,100.),(M, K, C))
  sol = solve(prob, DPRKN8(),tstops=0:1.0:100)

#These are the displacements at the free dof.
u3 = (x->x[7]).(sol.u)  #Node 1, θ
u4 = (x->x[8]).(sol.u)  #Node 2, x
u5 = (x->x[9]).(sol.u)  #Node 2, y
u6 = (x->x[10]).(sol.u) #Node 2, θ
u7 = (x->x[11]).(sol.u) #Node 3, x
u9 = (x->x[12]).(sol.u) #Node 3, θ

#Define the solution time vector.
t=sol.t

#Define some system stats.
num_timesteps = length(t)
num_nodes = size(node_geometry)[1]
num_dof = num_nodes * 3

#The displacements at the fixed dof are zero.
u1 = zeros(Float64, num_timesteps)
u2 = zeros(Float64, num_timesteps)
u8 = zeros(Float64, num_timesteps)

#Define the deformation field over time for plotting.

#Need [q1 q2 q3 q4] which is [y1 θ1 y2 θ2] for each beam element.

#Define global deformations for beam element #1.
q_e1 = Array{Array{Float64}}(undef, num_timesteps)

for i=1:num_timesteps

    q_e1[i] = [u2[i], u3[i], u5[i], u6[i]]

end

#Define global deformations for beam element #2.
q_e2 = Array{Array{Float64}}(undef, num_timesteps)

for i=1:num_timesteps

    q_e2[i] = [u5[i] u6[i] u8[i] u9[i]]

end

using Plots

#Define the beam shape function.
function calculate_element_shape(q, L, x, offset)
  
  a0 = q[1]
  a1 = q[2]
  a2 = 1/L^2*(-3*q[1]-2*q[2]*L+3*q[3]-q[4]*L)
  a3 = 1/L^3*(2*q[1]+q[2]*L-2*q[3]+q[4]*L)
  w = a0 .+ a1.*x .+ a2.*x.^2 .+ a3.*x.^3 .+ offset

  return w

end

#Calculate the element lengths
L = calculate_member_length(member_definitions, node_geometry)

@gif for i in 1:num_timesteps

  #Calculate the deformed shape of element 1.
  num_sub_elements = 20
  dx = L[1] / num_sub_elements
  x_e1 = 0.0:dx:L[1]
  offset = 0.0
  w_e1 = calculate_element_shape(q_e1[i], L[1], x_e1, offset)

  #Calculate the deformed shape of element 2.
  num_sub_elements = 20
  dx = L[2] / num_sub_elements
  x_e2 = 0:dx:(L[2])
  offset = 0.0
  w_e2 = calculate_element_shape(q_e2[i], L[2], x_e2, offset)

  #Define the beam shape.
  x = [x_e1[1:end]; maximum(x_e1) .+ x_e2[2:end]]
  beam_shape = [w_e1[1:end]; w_e2[2:end]]

  plot(x, beam_shape, legend=false, xlims=(0, 30), ylims=(-0.2, 0.2))

end






