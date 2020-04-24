using DifferentialEquations
using Dierckx  #for interpolating an aribtrary loading

#write equation of motion function
function skiddingcar!(du,u,p,t)

    m,pxSpline = p

    ptx=ptxSpline(t[1])   #define x direction load at time t

    #write the ODE m*d2ux/dx2=pxt as two first order ODEs

    #ux_dot=z        #velocity
    #z_dot=ptx/m     #acceleration

    #u[1] is ux
    #u[2] is z
    #du[1] is ux_dot
    #du[2] is z_dot

    du[1] = u[2]
    du[2] = ptx/m

end


#define mass of car
m = 1000


#define engine force
tx = 0:0.1:2.0 #seconds
a0=10000
ptx = a0.*ones(length(tx))
ptxSpline = Spline1D(tx,ptx)

uzero=[0., 0.]

tspan=(0.0,10.0)
prob = ODEProblem(skiddingcar!,uzero,tspan,(m,ptxSpline))
sol = solve(prob, tstops=0:0.1:10)

#extract solution
t=sol.t
ux_dot=(x->x[1]).(sol.u)
ux=(x->x[2]).(sol.u)

# using Plots
# plot(t, ux)


using Makie

CarWidth=3
CarLength=10

#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum(ux)-10.0
Ymin=Xmin
Xmax=maximum(ux)+10.0
Ymax=Xmax
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)

#define initial position
XYPos1=Node(ux[1])



#draw car
poly!(scene1, lift(xy->Point2f0[[CarLength/2+xy, CarWidth/2], [-CarLength/2+xy, CarWidth/2], [-CarLength/2+xy, -CarWidth/2], [CarLength/2+xy, -CarWidth/2]], XYPos1), color = :red, show_axis = false, limits=limits, overdraw =:true)


#draw ground and building centerline
# lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
# lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)


#animate

record(scene1, "animation.mp4", range(1, stop = length(ux), step=1)) do i

    if i==1
        sleep(1)
    end

    XYPos1[]=(ux[i])

    sleep(1/24)

end
