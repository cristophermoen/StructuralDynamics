#' # Lecture 4 - Simply supported beam vibration

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' Solve the PDE for simply supported beam vibration.

#' Define beam properties.
L = 100.0 #m
E = 1.99948e+11 #N/m^2
I = 385.0 * 25.4^4 /1000^4 #m^4, W14x38 
A = 11.2 * 25.4^2 / 1000^2 #m^2, W14x38
γ = 7700 #kg/m^3, density of steel

#' Define the mass of the beam per unit length.  
m = A * γ  #kg/m

#' Define the time range for the dynamic simulation.
t_start = 0.0 #sec
t_end = 10.0 #sec

#' Define the time discretization.
dt = 0.5 #sec

#' Define the beam discretization.
dx = 0.5  #m

#' We will solve $EI\frac{\partial^4 y}{\partial x^4} +  m\frac{\partial^2 y}{\partial t^2} = q(x,t)$ by using linear operators derived using finite differences.  The linear operators are developed 'by hand' for now.   The Julia language package (DiffEqOperators.jl)[https://github.com/SciML/DiffEqOperators.jl] has tools that can generate these operators more efficiently that what is done here.  Save that for later.

#' The PDE can be written as a linear equation $(EI[A] + m[B]){y} = {q}$ and solved as ${y} = (EI[A] + m[B])/{q}$.

#' Let's work on defining [A] first.

#' Calculate the number of segments, the number of nodes, and define the $x$ vector for the beam.
num_segments = floor(Int, L/dx)
num_nodes = num_segments + 1
x=0:L/num_segments:L 

#' Use the Fornberg algorithm to calculate derivative stencils.
include("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_4/helpers.jl")

#' Calculate the 4th derivative operator.  Use 5 nodes here for derivative accuracy. Calculate a centered difference here, at the middle of the stencil.
order = 4  
x_stencil = (0.0:1.0:4.0) .* dx  
x0_stencil = 2.0 .* dx 
stencil_4 = calculate_weights(order, x0_stencil, x_stencil)

#' For a simply-supported beam, the degrees of freedom for the first node and the last node are not considered in the solution because they are supports, i.e., y=0.  We don't need derivative stencils for these end nodes then.

#' We do need special stencils for the 2nd node and the (end - 1) nodes though, to impose the boundary condition $ \frac{d^2 y}{dx^2} = 0$ at the beam ends.  The approach proposed by [Souza](https://cimec.org.ar/ojs/index.php/mc/article/download/2662/2607) is used here.

bc_flag = 1 #simply supported
order = 4   #fourth derivative
stencil_4_left_boundary = calculate_boundary_stencils(bc_flag, dx, order)[1]
stencil_4_right_boundary = reverse(stencil_4_left_boundary)

#' Create the [A] operator.  We can define an [A] submatrix and then copy that through for each time step.  The [A] submatrix will have size num_nodes x num_nodes.
A_sub = zeros(Float64, (num_nodes, num_nodes))

for i = 1:num_nodes

    if i == 2

        A_sub[i,1:4] .= stencil_4_left_boundary

    elseif i == num_nodes - 1

        A_sub[i,end-3:end] .= stencil_4_right_boundary

    elseif (i != 1) & (i != num_nodes)

        A_sub[i,i-2:i+2] .= stencil_4

    end

end

#' The [A] operator has size of (num_nodes*num_timesteps) x (num_nodes*num_timesteps).
num_timesteps = floor(Int, (t_end-t_start)/dt + 1)
A = zeros(Float64, (num_nodes*num_timesteps, num_nodes*num_timesteps))

#' Copy A_sub through the [A] matrix for each time step.
for i=1:num_timesteps

    index = ((i-1)*num_nodes + 1):num_nodes*i

    A[index, index] .= A_sub

end

#' We need to remove the fixed (y=0) degrees of freedom from the solution. For this beam, that means nodes 1 and end.

#' Define the fixed dof in A_sub.
fixed_dof_sub = [1, num_nodes]

#' Define all the dof in A_sub.  
all_dof_sub = 1:num_nodes

#' Therefore the free dof in A_sub is...
free_dof_sub = setdiff(all_dof_sub, fixed_dof_sub)
num_free_dof_sub = length(free_dof_sub)

#' We need to remove the fixed dof at every timestep, so define that index vector.
free_dof = zeros(Int64, num_free_dof_sub*num_timesteps)

for i = 1:num_timesteps

    index = ((i-1)*num_free_dof_sub + 1):num_free_dof_sub*i
    free_dof[index] = free_dof_sub .+ (i-1)*num_nodes

end

#' [A_free] is the final form of the 4th order derivative operator [A] including boundary conditions.
A_free = A[free_dof, free_dof]

#' Now let's work on the [B] operator for the second derivative of y with respect to time t.


#' Calculate the 2nd derivative operator. Here 3 nodes are commonly used. Calculate a centered difference, at the middle of the stencil.
order = 2  
t_stencil = (0.0:1.0:2.0) .* dt  
t0_stencil = 1.0 .* dt 
stencil_2 = calculate_weights(order, t0_stencil, t_stencil)

#' We also need forward and backward facing stencils for when t=t_start and when t=t_end.
order = 2  
t_stencil = (0.0:1.0:2.0) .* dt  
t0_stencil = 0.0 .* dt 
stencil_2_forward = calculate_weights(order, t0_stencil, t_stencil)
stencil_2_backward = reverse(stencil_2_forward)

# #' Define B_sub
# B_sub = zeros(Float64, (num_timesteps, num_timesteps))

# for i = 1:num_timesteps

#     if i == 1

#         B_sub[i,1:3] .= stencil_2_forward

#     elseif i == num_timesteps

#         B_sub[i,end-2:end] .= stencil_2_backward

#     else

#         B_sub[i,i-1:i+1] .= stencil_2

#     end

# end


#' Define the operator matrix [B].
B = zeros(Float64, (num_nodes*num_timesteps, num_nodes*num_timesteps))

count = 1

for i = 1:num_timesteps

    for j = 1:num_nodes

        if i == 1
            column_index = (1:num_nodes:2*num_nodes+1) .+ (j-1)
            B[count, column_index] .= stencil_2_forward

        elseif i== num_timesteps
            column_index = (1:num_nodes:2*num_nodes+1) .+ (j-1) .+ num_nodes .* (i-3)
            B[count, column_index] .= stencil_2_backward

        else
            column_index = (1:num_nodes:2*num_nodes+1) .+ (j-1) .+ num_nodes .* (i-2)
            B[count, column_index] .= stencil_2

        end

        count = count + 1

    end

end

#' [B_free] is the final form of the 2th order derivative operator [B] including boundary conditions.
B_free = B[free_dof, free_dof]


#' Define the load on the beam $q(x,t)$.   
q = zeros(Float64, (num_nodes, num_timesteps))

#' Apply a uniform distributed load on the beam for the first second, then take it off.
q[:,1] .= 1.0  #N/m
q[:,2] .= 1.0  #N/m

#' Reshape $q(x,t) to prepare for the solution.
q = reshape(q, (num_nodes*num_timesteps, 1))

#' Use just the loads at the free dof in the solution.
q_free = q[free_dof]

#' Solve the linear system for the beam deflections over time.
y_free = (E .* I .* A_free .+ m .* B_free) \ q_free

#' Reshape the displacement vector to get a matrix of displacements, where the columns are the beam displacements shape over time.
y_free = reshape(y_free, (num_free_dof_sub, num_timesteps))

#' Add the zero displacements at the supports.
y =  [zeros(Float64, num_timesteps)'; y_free; zeros(Float64, num_timesteps)']



using Plots
plot(t_start:dt:t_end, y[26,:])