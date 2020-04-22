#Simulate the vibration of a simply supported beam with an initial imposed
#displacement at midspan.

cd(@__DIR__)  #change working directory to that of this file
include("MDOFDynamicsFunctions.jl")


                      #A Ix b(depth)
SectionProperties = [(16.0, 0.0, 50.0),   #Section property 1
                     (4.5, 0.0, 10.0),   #Section property 2
                     (43.98, 0.0, 400.0),   #Section property 3
                     (17.082, 0.0, 220)]   #Section property 4

#E  ρ (unit weight)    N/m^2  #kg/m^3 in this problem
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

#at DOF 4
t4 = 0:0.01:10 #seconds
a4 = 0.0
ω4 = 1.0
pt4 = a4*sin.(ω4*t4)
pt4Spline = Spline1D(t4,pt4)

#at DOF 5
t5 = 0:0.01:3 #seconds
a5 = 0.0
ω5 = 1.0
pt5 = a5*cos.(ω5*t5)
pt5Spline = Spline1D(t5,pt5)

#at DOF 6
t6 = 0:0.01:3 #seconds
a6 = 0.0
ω6 = 1.0
pt6 = a6*cos.(ω6*t6)
pt6Spline = Spline1D(t6,pt6)

#at DOF 7
t7 = 0:0.01:10 #seconds
a7 = 0.0
ω7 = 1.0
pt7 = a7*sin.(ω7*t7)
pt7Spline = Spline1D(t7,pt7)

#at DOF 8
t8 = 0:0.01:3 #seconds
a8 = 0.0
ω8 = 1.0
pt8 = a8*cos.(ω8*t8)
pt8Spline = Spline1D(t8,pt8)

#at DOF 9
t9 = 0:0.01:3 #seconds
a9 = 150.0
ω9 = 1.0
pt9 = a9*cos.(ω9*t)
pt9Spline = Spline1D(t9,pt9)

#at DOF 10
t10 = 0:0.01:3 #seconds
a10 = 0.0
ω10 = 1.0
pt10 = a10*cos.(ω10*t)
pt10Spline = Spline1D(t10,pt10)

#at DOF 11
t11 = 0:0.01:10 #seconds
a11 = 0.0
ω11 = 1.0
pt11 = a11*sin.(ω11*t11)
pt11Spline = Spline1D(t11,pt11)

#at DOF 12
t12 = 0:0.01:3 #seconds
a12 = 0.0
ω12 = 1.0
pt12 = a12*cos.(ω12*t12)
pt12Spline = Spline1D(t12,pt12)

#at DOF 13
t13 = 0:0.01:3 #seconds
a13 = 0.0
ω13 = 1.0
pt13 = a13*cos.(ω13*t13)
pt13Spline = Spline1D(t13,pt13)

#at DOF 14
t14 = 0:0.01:10 #seconds
a14 = 0.0
ω14 = 1.0
pt14 = a14*sin.(ω14*t14)
pt14Spline = Spline1D(t14,pt14)

#at DOF 15
t15 = 0:0.01:3 #seconds
a15 = 0.0
ω15 = 1.0
pt15 = a15*cos.(ω15*t15)
pt15Spline = Spline1D(t15,pt15)

#at DOF 16
t16 = 0:0.01:3 #seconds
a16 = 0.0
ω16 = 1.0
pt16 = a16*cos.(ω16*t)
pt16Spline = Spline1D(t16,pt16)


#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt2Spline, pt3Spline, pt4Spline, pt5Spline, pt6Spline, pt7Spline, pt8Spline, pt9Spline, pt10Spline, pt11Spline, pt12Spline, pt13Spline, pt14Spline, pt15Spline, pt16Spline = p

    pt1=pt1Spline(t[1])
    pt2=pt2Spline(t[1])
    pt3=pt3Spline(t[1])
    pt4=pt4Spline(t[1])
    pt5=pt5Spline(t[1])
    pt6=pt6Spline(t[1])
    pt7=pt7Spline(t[1])
    pt8=pt8Spline(t[1])
    pt9=pt9Spline(t[1])
    pt10=pt10Spline(t[1])
    pt11=pt11Spline(t[1])
    pt12=pt12Spline(t[1])
    pt13=pt13Spline(t[1])
    pt14=pt14Spline(t[1])
    pt15=pt15Spline(t[1])
    pt16=pt16Spline(t[1])
    pt=[pt1; pt2; pt3; pt4; pt5; pt6; pt7; pt8; pt9; pt10; pt11; pt12; pt13; pt14; pt15; pt16; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0]

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0                       u_0                         trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.; 0.], [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0], (0.,60.),(M, C, K, pt1Spline, pt2Spline, pt3Spline, pt4Spline, pt5Spline, pt6Spline, pt7Spline, pt8Spline, pt9Spline, pt10Spline, pt11Spline, pt12Spline, pt13Spline, pt14Spline, pt15Spline, pt16Spline))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:60)

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
u13_dot=(x->x[13]).(sol.u)
u14_dot=(x->x[14]).(sol.u)
u15_dot=(x->x[15]).(sol.u)
u16_dot=(x->x[16]).(sol.u)
u1=(x->x[17]).(sol.u)
u2=(x->x[18]).(sol.u)
u3=(x->x[19]).(sol.u)
u4=(x->x[20]).(sol.u)
u5=(x->x[21]).(sol.u)
u6=(x->x[22]).(sol.u)
u7=(x->x[23]).(sol.u)
u8=(x->x[24]).(sol.u)
u9=(x->x[25]).(sol.u)
u10=(x->x[26]).(sol.u)
u11=(x->x[27]).(sol.u)
u12=(x->x[28]).(sol.u)
u13=(x->x[29]).(sol.u)
u14=(x->x[30]).(sol.u)
u15=(x->x[31]).(sol.u)
u16=(x->x[32]).(sol.u)

t=sol.t


using Plots
plot(t,u9)
