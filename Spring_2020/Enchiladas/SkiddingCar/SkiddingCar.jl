using DifferentialEquations
using Dierckx  #for interpolating an aribtrary loading

#write equation of motion function
function skiddingcar!(du,u,p,t)

    m,Im,EngineForceSpline,SteeringSpline,c = p

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
    du[2] = pt*cos(θ)/m - c/m*du[1]

    #y direction
    du[3] = u[4]
    du[4] = pt*sin(θ)/m - c/m*du[3]

    #rotation
    du[5] = u[6]
    du[6] = mt/Im - c/m*du[5]

end

m1=1000   #mass of front of the car
m2=800  #mass of back of the car
r1=2.4   #distance from front car mass to center of gravity
r2=2.2   #distance from back car mass to center of gravity


c=0.00 #energy dissipation ratio

m = m1+m2   #total mass
Im= r1*m1+r2*m2   #rotational contribution of the mass

#define engine force
tEngineStart=0:0.01:1.0
EngineForceStart=zeros(length(tEngineStart))
tEngineDrive = 1.01:0.01:5.0 #seconds
a0=14500 #driving force, N
EngineForceDrive = a0.*ones(length(tEngineDrive))
tEngineEnd=5.01:0.01:10.0
EngineForceEnd=zeros(length(tEngineEnd))

tEngine=[tEngineStart;tEngineDrive;tEngineEnd]
EngineForce=[EngineForceStart;EngineForceDrive;EngineForceEnd]


plot(tEngine,EngineForce)

EngineForceSpline = Spline1D(tEngine,EngineForce)


plot(tEngine,EngineForceSpline.(tEngine))

#define steering
tNoSteering=0:0.01:3.5
tSteering=3.51:0.1:5.0
tAllSteering = [tNoSteering;tSteering]

a1=0.0 #constant for now,  sign change is to left or to the right
SteeringMoment = a1.*ones(length(tSteering))
AllSteeringMoment = [zeros(length(tNoSteering));SteeringMoment] #seconds

#create splines for interpolation
EngineForceSpline = Spline1D(tEngine,EngineForce)
SteeringSpline = Spline1D(tAllSteering,AllSteeringMoment)

#define initial conditions
         #ux0     #uy0     #uθ0
uzero=[0., 0., 0., 0., 0., 0.]

tspan=(0.0,5.0)
prob = ODEProblem(skiddingcar!,uzero,tspan,(m,Im,EngineForceSpline,SteeringSpline,c))
sol = solve(prob, tstops=0:0.01:5)

#extract solution
t=sol.t
ux=(x->x[1]).(sol.u)
ux_dot=(x->x[2]).(sol.u)
uy=(x->x[3]).(sol.u)
uy_dot=(x->x[4]).(sol.u)
uθ=(x->x[5]).(sol.u)
uθ=(x->x[6]).(sol.u)

ddux=diff(ux_dot)./diff(t)
ddux=[ddux[1];ddux]
dduy=diff(uy_dot)./diff(t)
dduy=[dduy[1];dduy]
dduθ=diff(uθ_dot)./diff(t)
dduθ=[dduθ[1];dduθ]


using Plots
plot(ux,uy, markershape =:circle)

plot(t, ux_dot)


plot(t,ux,xlabel="t[sec.]",ylabel="x-axis displacement[m]")
plot(t,uy,xlabel="t[sec.]",ylabel="y-axis displacement[m]")
plot(t,uθ,xlabel="t[sec.]",ylabel="rotation[rad]")

plot(t,ux_dot,xlabel="t[sec.]",ylabel="x-axis velocity[m/s]")
plot(t,uy_dot,xlabel="t[sec.]",ylabel="y-axis velocity[m/s]")
plot(t,uθ_dot,xlabel="t[sec.]",ylabel="rotation velocity[rad/s]")

plot(t,ddux,xlabel="t[sec.]",ylabel="x-axis accelartion[m/s^2]")
plot(t,dduy,xlabel="t[sec.]",ylabel="y-axis accelartion[m/s^2]")
plot(t,dduθ,xlabel="t[sec.]",ylabel="rotation accelartion[rad/s^2]")
