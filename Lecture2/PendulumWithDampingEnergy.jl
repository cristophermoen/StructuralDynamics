
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
dt=0.001 #seconds
totalt=10 #seconds

#define initial conditions
θo=pi/4
dθo=0.0

#solve the ODE
prob = SecondOrderODEProblem(PendulumEOM, [dθo], [θo], (0.,totalt), (g,L,c, m))
sol = solve(prob, DPRKN6(),tstops=0:dt:totalt)

#separate out θ(t),dθ(t), and t from solution
dθ=first.(sol.u)
θ=last.(sol.u)
t=sol.t


using Plots    #start Julia plot package


#calculate angular acceleration vs. time
    ddθ=diff(dθ)./diff(t)
    ddθ=[ddθ[1];ddθ]

#plot angular rotation, velocity, and acceleration
    p1=plot(t,θ,xaxis="t (sec.)",yaxis="\\theta")
    p2=plot(t,dθ,xaxis="t (sec.)",yaxis="d\\theta")
    p3=plot(t,ddθ,xaxis="t (sec.)",yaxis="dd\\theta")
    plot(p1,p2,p3,layout=(3,1),legend=false)

#let's study moment-rotation and energy

     #calculate moments acting first

        #calculate inertial moment
        momentI=m.*L.^2 .*ddθ
        #calculate gravity moment
        momentG=m.*g.*L*sin.(θ)
        #calculate damping moment
        momentD=c .*dθ*L.^2

     #plot M vs. time
        p1=plot(t,momentI,xaxis="t (sec.)",yaxis="momentI")
        p2=plot(t,momentG,xaxis="t (sec.)",yaxis="momentG")
        p3=plot(t,momentD,xaxis="t (sec.)",yaxis="momentD")
        plot(p1,p2,p3,layout=(3,1),legend=false)

    # plot M vs. θ
        p1=plot(θ,momentI,xaxis="rotation, \\theta (t) (radians)",yaxis="momentI, kN-m")
        p2=plot(θ,momentG,xaxis="rotation, \\theta (t) (radians)",yaxis="momentG, kN-m")
        p3=plot(θ,momentD,xaxis="rotation, \\theta (t) (radians)",yaxis="momentD, kN-m")
        plot(p1,p2,p3,layout=(3,1),legend=false)

    #calculate energy quantities
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
            p1=plot(t,energyI,xaxis="t (sec.)",yaxis="energyI, kN-m")
            p2=plot(t,energyG,xaxis="t (sec.)",yaxis="energyG, kN-m")
            p3=plot(t,energyD, xaxis="t (sec.)",yaxis="energyD, kN-m")
            plot(p1,p2,p3,layout=(3,1),legend=false)

    #calculate cumulative energy over time

            CumulativeEnergyI=accumulate(+, energyI)
            CumulativeEnergyG=accumulate(+, energyG)
            CumulativeEnergyD=accumulate(+, energyD)

            #add up all the energy
            #should be zero at any time t, however there is some numerical error
            TotalCEnergy=CumulativeEnergyI+CumulativeEnergyG+CumulativeEnergyD

            #plot cumulative energy vs. time
            p1=plot(t,CumulativeEnergyI,xaxis="t (sec.)",yaxis="CeI, kN-m")
            p2=plot(t,CumulativeEnergyG,xaxis="t (sec.)",yaxis="CeG, kN-m")
            p3=plot(t,CumulativeEnergyD,xaxis="t (sec.)",yaxis="CeD, kN-m")
            p4=plot(t,TotalCEnergy,xaxis="t (sec.)",yaxis="CeI+G+D, kN-m")
            plot(p1,p2,p3,p4,layout=(4,1),legend=false)
