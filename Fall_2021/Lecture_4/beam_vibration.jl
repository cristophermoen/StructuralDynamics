#' # Lecture 4 - Simply supported beam vibration

#' EN 560.630 Fall 2021

#' Department of Civil and Systems Engineering

#' Johns Hopkins University

#' Solve the PDE for simply supported beam vibration.

#' Define beam properties.
L = 5.0 #m
E = 1.99948e+11 #N/m^2
I = 1.0 #m^4
A = 0.1 #m^2
γ = 7700 #kg/m^3, density of steel

#' Define the mass of the beam per unit length.  
m = A * γ  #kg/m

#' Define the time range for the dynamic simulation.
t_start = 0.0 #sec
t_end = 10.0 #sec

#' Define the time discretization.
dt = 1.0 #sec

#' Define the beam discretization.
dx = 1.0  #m

#' We will solve $EI\frac{\partial^4 y}{\partial x^4} +  m\frac{\partial^2 y}{\partial t^2} = q(x,t)$ by using linear operators derived using finite differences.

#' The PDE can be written as a linear equation $(EI[A] + m[B]){y} = {q}$ and solved as ${y} = (EI[A] + m[B])/{q}$.

#' Let's work on defining [A] first.

#' Calculate the number of segments, the number of nodes, and define the $x$ vector for the beam.
num_segments = floor(Int, L/dx)
num_nodes = num_segments + 1
x=0:L/num_segments:L 

#' Use the Fornberg algorithm to calculate derivative stencils.
include("/Users/crismoen/.julia/dev/StructuralDynamics/Fall_2021/Lecture_4/helpers.jl")

#' Calculate the 4th derivative operator.  Use 5 nodes here to improve the derivative accuracy. Calculate a centered difference here, at the middle of the stencil.
order = 4  
x_stencil = (0.0:1.0:4.0) ./ dx  
x0_stencil = 2.0 ./ dx 
stencil_4 = calculate_weights(order, x0_stencil, x_stencil)

#' Calculate the 2nd derivative operator, 3 nodes are commonly used. Calculate a centered difference, at the middle of the stencil.
order = 2  
x_stencil = (0.0:1.0:2.0) ./ dx  
x0_stencil = 1.0 ./ dx 
stencil_2 = calculate_weights(order, x0_stencil, x_stencil)

#' Near the ends of the beam, a centered stencil won't work for the 4th derivative because we run out of nodes.  That's OK though, we can define a forward facing or backward facing stencil.
order = 4  
x_stencil = (0.0:1.0:4.0) ./ dx  
x0_stencil = 1.0 ./ dx 
stencil_4_forward_facing = calculate_weights(order, x0_stencil, x_stencil)
stencil_4_backward_facing = reverse(stencil_4_forward_facing)

#' Establish the beam boundary conditions.   We are assuming a simply supported beam and therefore $ \frac{d^2 y}{dx^2} = 0$ at each end.  We need a stencil that gives us this.

bc_flag = 1 #simply supported end
order = 2   #second derivative
stencil_left_end = calculate_boundary_stencils(bc_flag, dx, order)[1]
stencil_right_end = reverse(stencil_left_end)

#' Create the [A] operator.  Let's start with the portion of [A] at t=t_start.  This matrix will be num_nodes x num_nodes in size.
A_sub = zeros(Float64, (num_nodes, num_nodes))

for i = 1:num_nodes

    if i == 1

        A_sub[i,1:4] .= stencil_left_end

    elseif i == 2

        A_sub[i,1:5] .= stencil_4_forward_facing

    elseif i == num_nodes - 1

        A_sub[i,end-4:end] .= stencil_4_backward_facing

    elseif i == num_nodes

        A_sub[i,end-3:end] .= stencil_right_end

    else

        A_sub[i,i-2:i+2] .= stencil_4

    end

end
