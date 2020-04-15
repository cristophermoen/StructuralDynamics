
#Simulate the vibration of a simply supported beam with an initial imposed
#displacement at midspan.

cd(@__DIR__)  #change working directory to that of this file
include("MDOFDynamicsFunctions.jl")


                      #A Ix b(depth)
SectionProperties = [(1.0, 1.0, 2.0),   #Section property 1
                     (2.0, 2.4, 220)]   #Section property 2

#E  ρ (unit weight)    N/m^2  #kg/m^3 in this problem
MaterialProperties = [(1.0, 1.0),    #steel
                      (1E+5, 3000)]     #concrete

#x y
NodeGeometry = [0.0  0.0;
                12.0 0.0;
                24.0 0.0]

#iNode  jNode SectionProperties MaterialProperties
MemberDefinitions = [(1, 2, 1, 1),
                     (2, 3, 1, 1)]

          # node number, DOF u=1 v=2 ϕ=3
Supports = [(1,1),
            (1,2),
            (3,2)]



K = CalculateGlobalK(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

M = CalculateGlobalM(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

#calculate natural frequencies
#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
ωn_squared=eigvals(K,M)
ωn=sqrt.(ωn_squared)

#Raleigh damping
ωi=ωn[4]   #set first frequency anchor point i
ωj=ωn[5]   #set second frequency anchor point j
ζi=0.05    #modal viscous damping ratio in mode i
ζj=0.05    #modal viscous damping ratio in mode j

a0,a1=2*inv([1/ωi ωi;1/ωj ωj])*[ζi;ζj]  #Eq. 11.4.9 from Chopra
C=a0*M+a1*K


#define arbitrary forcing functions
#at DOF 1
t1 = 0:0.01:10 #seconds
a1 = 0.0
ω1 = 1.0
pt1 = a1*sin.(ω1*t1)
pt1Spline = Spline1D(t1,pt1)

#at DOF 2
t2 = 0:0.01:3 #seconds
a2 = 0.0
ω2 = 1.0
pt2 = a2*cos.(ω2*t2)
pt2Spline = Spline1D(t2,pt2)

#at DOF 3
t3 = 0:0.01:3 #seconds
a3 = 0.0
ω3 = 1.0
pt3 = a2*cos.(ω3*t3)
pt3Spline = Spline1D(t3,pt3)

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt2Spline, pt3Spline = p

    pt1=pt1Spline(t[1])
    pt2=pt2Spline(t[1])
    pt3=pt3Spline(t[1])
    pt=[pt1; pt2; pt3; 0.0; 0.0; 0.0]

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0                       u_0                         trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.; 0.; 0.; 0.], [0.0; 0.0; 0.1; 0.0; 0.0; 0.0], (0.,60.),(M, C, K, pt1Spline, pt2Spline, pt3Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:60)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u3_dot=(x->x[3]).(sol.u)
u4_dot=(x->x[4]).(sol.u)
u5_dot=(x->x[5]).(sol.u)
u6_dot=(x->x[6]).(sol.u)
u1=(x->x[7]).(sol.u)
u2=(x->x[8]).(sol.u)
u3=(x->x[9]).(sol.u)
u4=(x->x[10]).(sol.u)
u5=(x->x[11]).(sol.u)
u6=(x->x[12]).(sol.u)

t=sol.t


using Plots
plot(t,u3)
