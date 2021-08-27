#Simulate the vibration of a cantilever beam with a harmonic force

cd(@__DIR__)  #change working directory to that of this file
include("EnchiladaFunctions.jl")

                       #Area   #Ix      #b
SectionProperties = [(0.023, 4.76e-6, 0.457)]   #50/50 Aluminum and Fiberglass

#E  ρ (unit weight)    N/m^2     kg/m^3
MaterialProperties = [(7.065e10, 2179)]    #50/50 Aluminum and Fiberglass

                  #x y
NodeGeometry = [0.0 0.0;
                0.458 0.0;
                0.915 0.0;
                1.37 0.0;
                1.83 0.0]

#iNode  jNode SectionProperties MaterialProperties
MemberDefinitions = [(1, 2, 1, 1),
                     (2, 3, 1, 1),
                     (3, 4, 1, 1),
                     (4, 5, 1, 1)]

          # node number, DOF u=1 v=2 ϕ=3
Supports = [(1,1),
            (1,2),
            (1,3)]

K, FreeDOF = CalculateGlobalK(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

M, FreeDOF = CalculateGlobalM(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

#calculate natural frequencies
#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
ωn_squared=eigvals(K,M)
ωn=sqrt.(ωn_squared)

#Raleigh damping
ωi=ωn[1]   #set first frequency anchor point i
ωj=ωn[5]   #set second frequency anchor point j
ζi=0.02    #modal viscous damping ratio in mode i
ζj=0.02    #modal viscous damping ratio in mode j

a0,a1=2*inv([1/ωi ωi;1/ωj ωj])*[ζi;ζj]  #Eq. 11.4.9 from Chopra
C=a0*M+a1*K


#Define person load at FDOF 11 on node 5
f1 = 0.5 #Hertz
t1 = 0:0.01:0.5 #Seconds
a1 = -5422   #Newtons
ω1 = 2*pi*f1
tFirst = 0:0.01:0.01 #Seconds
ptFirst = zeros(length(tFirst)) #First Force Zero
tSecond = 0.02:0.01:0.52 #Seconds
ptSecond = a1*sin.(ω1*tSecond) #Second Force Applied
tThird = 0.53:0.01:10.0 #Seconds
ptThird=zeros(length(tThird)) #Third Force Zero
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt1Spline = Spline1D(tAll,ptAll)
pt1DOF=findfirst(x->x==11, FreeDOF)  #Apply Load at FDOF 11

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt1DOF = p

    pt=zeros(size(M)[1])  #define vector of all the forcing functions at each free DOF

    pt[pt1DOF]=pt1Spline(t[1])  #assign load to node 5, vertical DOF
    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0             u_0                                                            trange
  prob = SecondOrderODEProblem(mdof, zeros(size(M)[1]), zeros(size(M)[1]), (0.,10.),(M, C, K, pt1Spline, pt1DOF))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

  t=sol.t

  u1_dot=(x->x[1]).(sol.u)
  u2_dot=(x->x[2]).(sol.u)
  u3_dot=(x->x[3]).(sol.u)
  u4_dot=(x->x[4]).(sol.u)
  u5_dot=(x->x[5]).(sol.u)
  u6_dot=(x->x[6]).(sol.u)
  u7_dot=(x->x[7]).(sol.u)
  u8_dot=(x->x[8]).(sol.u)
  u9_dot=(x->x[9]).(sol.u)
  u10_dot=(x->x[10]).(sol.u)
  u11_dot=(x->x[11]).(sol.u)
  u12_dot=(x->x[12]).(sol.u)

  u1=(x->x[13]).(sol.u)
  u2=(x->x[14]).(sol.u)
  u3=(x->x[15]).(sol.u)
  u4=(x->x[16]).(sol.u)
  u5=(x->x[17]).(sol.u)
  u6=(x->x[18]).(sol.u)
  u7=(x->x[19]).(sol.u)
  u8=(x->x[20]).(sol.u)
  u9=(x->x[21]).(sol.u)
  u10=(x->x[22]).(sol.u)
  u11=(x->x[23]).(sol.u)
  u12=(x->x[24]).(sol.u)


using Plots
plot(t,u11,xlabel="t[sec]",ylabel="displacement[meters]",title="Displacement vs. Time at End of Diving Board")  #plot displacement at DOF 26
