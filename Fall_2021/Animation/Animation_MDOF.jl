using LaTeXStrings, ProgressMeter, Interpolations

function BeamShape(q1,q2,q3,q4,L,x,offset)
    a0=q1
    a1=q2
    a2=1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3=1/L^3*(2*q1+q2*L-2*q3+q4*L)
    w=a0 .+ a1.*x .+ a2.*x.^2 .+ a3.*x.^3 .+ offset
end

function anim_plot(t_eq, displacements, perm_deformations, ground_motion, h_floors, DOF_disp_plot, anim_settings, bldg_representation = "column", d_shape = 20)
    # t_eq == time array inputed for the solver
    # displacements == (t, u)
    # perm_deformations == (perm_t, perm_u) or () if they are not included
    # ground_motion == (ground_t, ground_u) or () if they are not included
    
        # t == Array of time from the solver
        # u == [[u_dof1_array], [u_dof2_array], ..., [u_dofn_array]] -> Array of Arrays of different displacements for each DOF (if only 1 DOF, input should be [u])
        # perm_t == Array of time that aligns with the perminant deformations
        # perm_u == [[perm_u_dof1_array], [perm_u_dof2_array], ..., [perm_u_dofn_array]] -> Array of Arrays of different displacements for each DOF (if only 1 DOF, input should be [perm_u])
        # ground_t == Array of time that aligns with the ground motion
        # ground_u == Array of ground motion displacements

    # h_floors == [h_1, h_2, ..., h_n] -> Array of floor heights (if only 1 floor, input should be [h_1])
    # DOF_disp_plot == Array of the DOFs to be tracked in the second plot (use 0 for ground motion)
    # displacements == (fps, file_name) or () for default
    # bldg_representation == "column" or "frame" if you want to see building represented as a single column or as a frame (leave blank for column as default)
    # d_shape == Number of points between each floor for the deformation plot (leave blank for default)
                
        # fps == Number frames per second -> (recomend 10 to 20)
        # file_name == "file_name.gif" or "file_name.mp4"

    t, u = displacements

    if length(anim_settings) > 0
        fps, file_name, speed = anim_settings
    else
        fps, file_name, speed = (20, "Animate.mp4",1)
    end
    
    if length(perm_deformations) > 0
        perm_t, perm_u = perm_deformations
        perm_u_merged = Vector{Vector{Float64}}(undef, length(perm_u[1]))
        for j in eachindex(perm_u[1])
            perm_u_merged[j] = [perm_u[1][j]]
            if length(perm_u) > 1
                for i in eachindex(perm_u[2:end])
                    perm_u_merged[j] = [perm_u_merged[j];perm_u[i+1][j]]
                end
            end
        end
        
        t_index1 = zeros(length(perm_t))
        for ike in eachindex1(perm_t)
            t_index1[ike] = ike
        end
        time2index1 = Spline1D(perm_t,t_index1)
        t_index2u1 = interpolate(perm_u_merged, BSpline(Linear()))
        perm_u_dict(t_anim) = t_index2u1[time2index1(t_anim)]
    end
    
    if length(ground_motion) > 0
        ground_t, ground_u = ground_motion
    else
        ground_t = t
        ground_u = zeros(length(u[1]))
    end
    ground_u_dict = Spline1D(ground_t, ground_u)
    
    t_anim = (t_eq[1]:1/fps:t_eq[end])
    u_merged = Vector{Vector{Float64}}(undef, length(u[1]))
    for j in eachindex(u[1])
        u_merged[j] = [u[1][j]]
        if length(u) > 1
            for i in eachindex(u[2:end])
                u_merged[j] = [u_merged[j];u[i+1][j]]
            end
        end
    end
    
    t_index2 = zeros(length(t))
    for ike in eachindex2(t)
        t_index2[ike] = ike
    end
    time2index2 = Spline1D(t,t_index2)
    t_index2u2 = interpolate(u_merged, BSpline(Linear()))
    u_dict(t_anim) = t_index2u2[time2index2(t_anim)]

    
    color_a = colorant"rgb(8,156,252)"
    color_b = colorant"rgb(232,108,68)"
    color_c = :grey35
    
    x_max = 0
    x_min = 0

    for i in eachindex(t_anim)
        x_max = max(x_max, maximum(ground_u_dict(t_anim[i]) .+ [[0];u_dict(t_anim[i])]))
        x_min = min(x_min, minimum(ground_u_dict(t_anim[i]) .+ [[0];u_dict(t_anim[i])]))
    end
    
    if bldg_representation == "frame" || bldg_representation == "Frame"
        col = [-(x_max - x_min)*0.2,(x_max - x_min)*0.2]
    elseif bldg_representation == "column" || bldg_representation == "Column"
        col = [0]
    end
    
    t_anim_end = eachindex(t_anim)[end]
    prog = Progress(t_anim_end, "Building Animation: ")
    
    anim = @animate for i in (1:t_anim_end)
        ys = [[0];h_floors]
        xs = ground_u_dict(t_anim[i]) .+ [[0];u_dict(t_anim[i])]
        
        xvals = []
        yvals = []

        x1 = Vector{Vector{Float64}}(undef,length(xs))
        y1 = Vector{Vector{Float64}}(undef,length(ys))
        
        if length(perm_deformations) > 0
            perm_xs = ground_u_dict(t_anim[i]) .+ [[0];perm_u_dict(t_anim[i])]
            perm_xvals = []
            perm_x1 = Vector{Vector{Float64}}(undef,length(perm_xs))
        end
        
        for n in eachindex(ys[1:end-1])
            x1[n] = zeros(d_shape+1)
            y1[n] = zeros(d_shape+1)
            
            h = ys[n+1] - ys[n]
            Δ = xs[n+1] - xs[n]
            
            if length(perm_deformations) > 0
                perm_x1[n] = zeros(d_shape+1)
                perm_Δ = perm_xs[n+1] - perm_xs[n]
            end
            
            for m in (0:d_shape)
                x1[n][m+1] = BeamShape(0,0,1,0,h,h*m/d_shape,0) * Δ
                y1[n][m+1] = h*m/d_shape

                if length(perm_deformations) > 0
                    perm_x1[n][m+1] = BeamShape(0,0,1,0,h,h*m/d_shape,0) * perm_Δ
                end
            end    
            xvals = [xvals;x1[n] .+ (xs[n])]
            yvals = [yvals;y1[n] .+ (ys[n])]
            
            if length(perm_deformations) > 0
                perm_xvals = [perm_xvals;perm_x1[n] .+ (perm_xs[n])]
            end
        end
        
        plt1 = plot(xlabel=L"\mathrm{Displacement \quad [m]}", ylabel=L"\mathrm{Height \quad [m]}", 
            legend=false, xlims = ((x_min + col[1])*1.1,(x_max+ col[end])*1.1),ylims = (0,h_floors[end]*1.02), title = "\$t = $(round(t_anim[i], digits=1))\$")
        if length(perm_deformations) > 0
            for a in eachindex(col)
                plt1 = plot!(perm_xvals .+ col[a],yvals, linestyle = :dash, color = color_a)
            end
            if length(col) > 1
                for a in eachindex(h_floors)
                    plt1 = plot!([col[1],col[end]] .+ ground_u_dict(t_anim[i]) .+ perm_u_dict(t_anim[i])[a],[h_floors[a],h_floors[a]], color = color_a,linestyle = :dash)
                end
            end
            plt1 = scatter!([ground_u_dict(t_anim[i]) .+ perm_u_dict(t_anim[i])],h_floors, markerstrokewidth=0, color = color_a)
        end
        
        for a in eachindex(col)
            plt1 = plot!(xvals .+ col[a],yvals, color = color_b)
        end
        if length(col) > 1
            for a in eachindex(h_floors)
                plt1 = plot!([col[1],col[end]] .+ ground_u_dict(t_anim[i]) .+ u_dict(t_anim[i])[a],[h_floors[a],h_floors[a]], color = color_b)
            end
        end
        plt1 = scatter!(ground_u_dict(t_anim[i]) .+ u_dict(t_anim[i]),h_floors, color = color_b)
        if length(col) > 1 && length(ground_motion) > 0
            plt1 = plot!([ground_u_dict(t_anim[i]),ground_u_dict(t_anim[i])],[h_floors[end]/100,0],arrow=(:closed, 1),color=:black,linewidth=0.5,label="")
        end
        
        plt2_sub = Array{Any}(nothing, length(DOF_disp_plot))
        for ike in eachindex(DOF_disp_plot)
            plt2_sub[ike] = plot( ylabel=L"\mathrm{u \quad [m]}") #xlabel=L"\mathrm{Time \quad [sec.]}",
            if DOF_disp_plot[ike] == 0
                plt2_sub[ike] = plot!(t,ground_u_dict(t), label = L"\textrm{Ground\; Motion}", color = color_b)
                plt2_sub[ike] = plot!([t_anim[i],t_anim[i]],[minimum(ground_u)*1.1,maximum(ground_u)*1.1],label = false,linestyle = :dash, color = color_c)
                plt2_sub[ike] = scatter!([t_anim[i]],[ground_u_dict(t_anim[i])],label = false , markerstrokewidth=0, color = color_c)   
            else
                if length(perm_deformations) > 0
                    plt2_sub[ike] = plot!(perm_t,perm_u[DOF_disp_plot[ike]], label = "Perminant Deformation", color = color_a)
                end
                plt2_sub[ike] = plot!(t,u[DOF_disp_plot[ike]], label = L"\textrm{DOF\; %$(DOF_disp_plot[ike])}", color = color_b)
                plt2_sub[ike] = plot!([t_anim[i],t_anim[i]],[minimum(u[DOF_disp_plot[ike]])*1.1,maximum(u[DOF_disp_plot[ike]])*1.1],linestyle = :dash, label = false, color = color_c)
                plt2_sub[ike] = scatter!([t_anim[i]],[u_dict(t_anim[i])[DOF_disp_plot[ike]]],label = false , markerstrokewidth=0, color = color_c)   
            end
        end
        plt2_sub[end] = plot!(xlabel=L"\mathrm{Time \quad [sec.]}",bottom_margin = 20Plots.mm)
        plt2 = plot(plt2_sub..., layout=(length(DOF_disp_plot),1))
        if length(DOF_disp_plot) <= 3
            plot(plt1,plt2,layout = grid(1, 2, widths=[0.4 ,0.6]),size = (1200,450),left_margin = 20Plots.mm)
        else
            plt1_sized = plot(plt1,plot(legend=false,grid=false,foreground_color_subplot=:white), layout = grid(2, 1, heights=[3/length(DOF_disp_plot), 1-3/length(DOF_disp_plot)]))
            plot(plt1_sized,plt2,layout = grid(1, 2, widths=[0.4 ,0.6]),size = (1200,200*length(DOF_disp_plot)),left_margin = 20Plots.mm)
        end
        next!(prog)
    end
    gif(anim, file_name, fps = fps*speed)
end
