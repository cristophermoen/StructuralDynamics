#Simulate the vibration of a simply supported beam with an initial imposed
#displacement at midspan.

cd(@__DIR__)  #change working directory to that of this file
include("MDOFDynamicsFunctions.jl")


                      #A Ix b(depth)
SectionProperties = [(16.0, 0.0, 50.0),   #Section property 1
                     (4.5, 0.0, 10.0),   #Section property 2
                     (43.98, 0.0, 400.0),   #Section property 3
                     (17.082, 0.0, 220)]   #Section property 4

#E  ρ (unit weight)    ksi      lb/ft^3 in this problem
MaterialProperties = [(28000000, 485.0),    #Wrought Iron
                      (1.5E+7, 450.0)]     #Cast Iron

#x y
NodeGeometry = [0.0  0.0;   #Node 1
                138.0 154.0;    #Node 2
                138.0 0.0;   #Node 3
                276.0 154.0;   #Node 4
                276.0 0.0;   #Node 5
                414.0 154.0;   #Node 6
                414.0 0.0;    #Node 7
                552.0 154.0;   #Node 8
                552.0 0.0;   #Node 9
                690.0 154.0;   #Node 10
                690.0 0.0;   #Node 11
                828.0 154.0;   #Node 12
                828.0 0.0;   #Node 13
                966.0 154.0;   #Node 14
                966.0 0.0;   #Node 15
                1104.0 0.0]   #Node 16

#iNode  jNode SectionProperties MaterialProperties
MemberDefinitions = [(1, 2, 3, 2),      #Mbr 1
                     (1, 3, 1, 1),      #Mbr 2
                     (3, 2, 4, 2),      #Mbr 3
                     (2, 4, 3, 2),      #Mbr 4
                     (2, 5, 2, 1),      #Mbr 5
                     (3, 5, 1, 1),      #Mbr 6
                     (5, 4, 4, 2),      #Mbr 7
                     (4, 6, 3, 2),      #Mbr 8
                     (4, 7, 2, 1),      #Mbr 9
                     (5, 6, 2, 1),      #Mbr 10
                     (5, 7, 1, 1),      #Mbr 11
                     (7, 6, 4, 2),      #Mbr 12
                     (6, 8, 3, 2),      #Mbr 13
                     (6, 9, 2, 1),      #Mbr 14
                     (7, 8, 2, 1),      #Mbr 15
                     (7, 9, 1, 1),      #Mbr 16
                     (9, 8, 4, 2),      #Mbr 17
                     (8, 10, 3, 2),      #Mbr 18
                     (8, 11, 2, 1),      #Mbr 19
                     (9, 10, 2, 1),      #Mbr 20
                     (9, 11, 1, 1),      #Mbr 21
                     (11, 10, 4, 2),      #Mbr 22
                     (10, 12, 3, 2),      #Mbr 23
                     (10, 13, 2, 1),      #Mbr 24
                     (11, 12, 2, 1),      #Mbr 25
                     (11, 13, 1, 1),      #Mbr 26
                     (13, 12, 3, 2),      #Mbr 27
                     (12, 14, 3, 2),      #Mbr 28
                     (13, 14, 2, 1),      #Mbr 29
                     (13, 15, 1, 1),      #Mbr 30
                     (15, 14, 3, 2),      #Mbr 31
                     (14, 16, 3, 2),      #Mbr 32
                     (15, 16, 1, 1)]      #Mbr 33

          # node number, DOF u=1 v=2 ϕ=3
Supports = [(1,1),
            (1,2),
            (1,3),
            (2,3),
            (3,3),
            (4,3),
            (5,3),
            (6,3),
            (7,3),
            (8,3),
            (9,3),
            (10,3),
            (11,3),
            (12,3),
            (13,3),
            (14,3),
            (15,3),
            (16,3),
            (16,2)]


K, FreeDOF = CalculateGlobalK(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

M, FreeDOF = CalculateGlobalM(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

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


#define loads

#define pedestrian load at node 9
#the vertical DOF number for node 9 is 26
#at DOF 26
t1 = 0:0.01:10 #seconds
a1 = 1.0
ω1 = 1.0
pt1 = a1*sin.(ω1*t1)
pt1Spline = Spline1D(t1,pt1)
pt1DOF=findfirst(x->x==26, FreeDOF)  #vector location of load

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline,pt1DOF = p

    pt=zeros(size(M)[1])  #define vector of all the pts at each free DOF

    pt[pt1DOF]=pt1Spline(t[1])  #assign load to node 9, vertical DOF

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0       u_0         trange
  prob = SecondOrderODEProblem(mdof, zeros(size(M)[1]), zeros(size(M)[1]), (0.,60.),(M, C, K, pt1Spline, pt1DOF))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:60)

  t=sol.t

TotalNumberOfDOF=size(NodeGeometry)[1]*3
u_dot=zeros(length(t),TotalNumberOfDOF)
u=zeros(length(t),TotalNumberOfDOF)

for i=1:length(FreeDOF)
    u_dot[:,FreeDOF[i]]=(x->x[i]).(sol.u)
end

for i=1:length(FreeDOF)
    Increment=length(FreeDOF)
    u[:,FreeDOF[i]]=(x->x[i+Increment]).(sol.u)
end


using Plots
plot(t,u[:,26])  #plot displacement at DOF 26
