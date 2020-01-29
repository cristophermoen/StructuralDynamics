
#define Julia ODE package that we will use as solver
using OrdinaryDiffEq

#define constants
m=0.035    #kg
g=9.8  #m/s^2
L=0.5 #meters
c=0.018  #wind resistance parameter

#write equation of motion for a pendulum
function PendulumEOM(ddθ,dθ,θ,p,t)
    g, L, c, m = p
    ddθ[1] = -g*sin(θ[1])/L - c*dθ[1]/m
end

#solver controls
dt=0.01 #seconds
totalt=20 #seconds

#define initial conditions
θo=pi/4
dθo=0.0

#solve the ODE
prob = SecondOrderODEProblem(PendulumEOM, [θo], [dθo], (0.,totalt), (g,L,c, m))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out θ(t),dθ(t), and t from solution
θ=first.(sol.u)
dθ=last.(sol.u)   #angular velocity
t=sol.t

# #plot results
using Plots    #start Julia plot package
# plot(t,θ,linewidth=1,title="pendulum free vibration",
#     xaxis="time, t (sec.)",yaxis="rotation, \\theta (t) (radians)", legend=false)

#let's study load deformation and energy

    #calculate forces acting first

        #calculate intertial force
        #need acceleration, take numerical derivative of velocity
        ddθ=[0; diff(dθ)]./[0; diff(t)]
        fI=m.*L.^2 .*ddθ
        fI[1]=0   #clean up numerical messiness
        #calculate gravity force
        fG=m.*g.*L*sin.(θ)
        #calculate damping force
        fD=c .*dθ

    #now calculate moments

        #inertial moment
        momentI=fI.*L
        #gravity moment
        momentG=fG.*L
        #damping moment
        momentD=fD.*L

    #plot M vs. θ
        p1=plot(θ,momentI,xaxis="rotation, \\theta (t) (radians)",yaxis="momentI, kN-m")
        p2=plot(θ,momentG,xaxis="rotation, \\theta (t) (radians)",yaxis="momentG, kN-m")
        p3=plot(θ,momentD,xaxis="rotation, \\theta (t) (radians)",yaxis="momentD, kN-m")
        plot(p1,p2,p3,layout=(3,1),legend=false)

    #calculate energy quantities
            # using Statistics

            #write a little function to help with this
            #integrate little chunks of M vs. θ
            function CalculateEnergy(M,θ)
                energy=zeros(length(M))
                for i=1:length(M)
                    energy[i]=M[i]*θ[i]
                end
                return energy
            end

            #inertial energy
            energyI=CalculateEnergy(momentI,θ)

            #potential energy of mass
            energyG=CalculateEnergy(momentG,θ)

            # energy dissipated by air resistance
            energyD=CalculateEnergy(momentD,θ)

            #plot instantaneous energy vs. time
            p1=plot(t,energyI)
            p2=plot(t,energyG)
            p3=plot(t,energyD)

    #calculate cumulative energy

            CumulativeEnergyI=accumulate(+, energyI)
            CumulativeEnergyG=accumulate(+, energyG)
            CumulativeEnergyD=accumulate(+, energyD)

            #plot cumulative energy vs. time
            p1=plot(t,CumulativeEnergyI)
            p2=plot(t,CumulativeEnergyG)
            p3=plot(t,CumulativeEnergyD)
            plot(p1,p2,p3,layout=(3,1),legend=false)
