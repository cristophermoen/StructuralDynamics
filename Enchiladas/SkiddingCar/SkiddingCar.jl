using DifferentialEquations
using Dierckx  #for interpolating an aribtrary loading

#write equation of motion function
function skiddingcar!(du,u,p,t)

    m,Im,EngineForceSpline,SteeringSpline = p

    pt=EngineForceSpline(t[1])   #define x direction load at time t

    mt=SteeringSpline(t[1])

    #write the ODE m*d2ux/dx2=pxt as two first order ODEs

    #ux_dot=z        #velocity
    #z_dot=ptx/m     #acceleration

    #u[1] is ux
    #u[2] is z
    #du[1] is ux_dot
    #du[2] is z_dot

    θ=u[6]

    #x direction
    du[1] = u[2]
    du[2] = pt*cos(θ)/m

    #y direction
    du[3] = u[4]
    du[4] = pt*sin(θ)/m

    #rotation
    du[5] = u[6]
    du[6] = mt/Im

end

m1=500   #mass of front of the car
m2=500   #mass of back of the car
r1=1.2   #distance from front car mass to center of gravity
r2=1.2   #distance from back car mass to center of gravity

m = m1+m2   #total mass
Im= r1*m1+r2*m2   #rotational contribution of the mass

#define engine force
tEngine = 0:0.1:2.0 #seconds
a0=10000  #constant for now
EngineForce = a0.*ones(length(tEngine))
EngineForceSpline = Spline1D(tEngine,EngineForce)

#define steering
tNoSteering=0:0.1:2.0
tSteering=2.1:0.1:3.0
tAllSteering = [tNoSteering;tSteering]

a1=100  #constant for now,  sign change is to left or to the right
SteeringMoment = a1.*ones(length(tSteering))
AllSteeringMoment = [zeros(length(tNoSteering));SteeringMoment] #seconds

#create splines for interpolation
EngineForceSpline = Spline1D(tEngine,EngineForce)
SteeringSpline = Spline1D(tAllSteering,AllSteeringMoment)

#define initial conditions
uzero=[0., 0., 0., 0., 0., 0.]

tspan=(0.0,10.0)
prob = ODEProblem(skiddingcar!,uzero,tspan,(m,Im,EngineForceSpline,SteeringSpline))
sol = solve(prob, tstops=0:0.1:10)

#extract solution
t=sol.t
ux_dot=(x->x[1]).(sol.u)
ux=(x->x[2]).(sol.u)
uy_dot=(x->x[3]).(sol.u)
uy=(x->x[4]).(sol.u)
uθ_dot=(x->x[5]).(sol.u)
uθ=(x->x[6]).(sol.u)


using Plots
plot(ux,uy, markershape =:circle)




# using Makie
#
# CarWidth=3
# CarLength=10
#
# #define window and limits
# scene1=Scene(resolution=(1000,1000))
# Xmin=minimum(ux)-10.0
# Ymin=Xmin
# Xmax=maximum(ux)+10.0
# Ymax=Xmax
# Xrange=Xmax-Xmin
# Yrange=Ymax-Ymin
# limits = FRect(Xmin, Ymin, Xrange, Yrange)
#
# #define initial position
# XYPos1=Node(ux[1])
#
#
#
# #draw car
# poly!(scene1, lift(xy->Point2f0[[CarLength/2+xy, CarWidth/2], [-CarLength/2+xy, CarWidth/2], [-CarLength/2+xy, -CarWidth/2], [CarLength/2+xy, -CarWidth/2]], XYPos1), color = :red, show_axis = false, limits=limits, overdraw =:true)
#
#
# #draw ground and building centerline
# # lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
# # lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)
#
#
# #animate
#
# record(scene1, "animation.mp4", range(1, stop = length(ux), step=1)) do i
#
#     if i==1
#         sleep(1)
#     end
#
#     XYPos1[]=(ux[i])
#
#     sleep(1/24)
#
# end
