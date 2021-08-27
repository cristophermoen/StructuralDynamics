#Simulate the vibration of a simply supported beam with an initial imposed
#displacement at midspan.

cd(@__DIR__)  #change working directory to that of this file
include("MDOFDynamicsFunctions.jl")


                  #A(in^2) Ix b(depth)
SectionProperties = [(16.0, 0.0, 50.0),   #Section property 1
                     (4.5, 0.0, 10.0),   #Section property 2
                     (43.98, 0.0, 400.0),   #Section property 3
                     (17.082, 0.0, 220)]   #Section property 4

#E  ρ (unit weight)    lb/in^2      lb/in^3    in this problem
MaterialProperties = [(28000000, 485.0*0.000578704),    #Wrought Iron
                      (1.5E+7, 450.0*0.000578704)]     #Cast Iron

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
                     (3, 2, 2, 1),      #Mbr 3
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
                     (15, 14, 2, 1),      #Mbr 31
                     (14, 16, 3, 2),      #Mbr 32
                     (15, 16, 1, 1)]      #Mbr 33

          # node number, DOF u=1 v=2 ϕ=3
Supports = [(1,1),   #Pin Support
            (1,2),  #Pin Support
            (1,3),  #No moment
            (2,3),  #No moment
            (3,3),  #No moment
            (4,3),  #No moment
            (5,3),  #No moment
            (6,3),  #No moment
            (7,3),  #No moment
            (8,3),  #No moment
            (9,3),  #No moment
            (10,3), #No moment
            (11,3), #No moment
            (12,3), #No moment
            (13,3), #No moment
            (14,3), #No moment
            (15,3), #No moment
            (16,3), #No moment
            (16,2)] #Roller Support


K, FreeDOF = CalculateGlobalK(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

M, FreeDOF = CalculateGlobalM(MemberDefinitions, SectionProperties, MaterialProperties, NodeGeometry, Supports)

#calculate natural frequencies
#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
ωn_squared=eigvals(K,M)
ωn=sqrt.(ωn_squared)

#Raleigh damping
ωi=ωn[4]   #set first frequency anchor point i
ωj=ωn[5]   #set second frequency anchor point j
ζi=0.01    #modal viscous damping ratio in mode i
ζj=0.01    #modal viscous damping ratio in mode j

a0,a1=2*inv([1/ωi ωi;1/ωj ωj])*[ζi;ζj]  #Eq. 11.4.9 from Chopra
C=a0*M+a1*K


#define loads
f1 = 2.3 #hz  Equals the rate of marching of a light infantry unit is 140 paces/minute
a1 = 3*6*(170+60)   #3 rows of soldiers, 6 column deep with an average weight of 170l plus combat load
ω1 = 2*pi*f1        #Natural frequency of the force applied


#define pedestrian load at node 3
#the vertical DOF number for node 3 is 8
#Apply load starting from zero to 70 seconds
tFirst = 0:0.01:0.01 #seconds
ptFirst=zeros(length(tFirst))
tSecond=0.02:0.01:70 #seconds  Turn on the load
ptSecond=a1*sin.(ω1*tSecond)
tThird=70.01:0.01:130          #Turn off the load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt2Spline = Spline1D(tAll,ptAll)
pt2DOF=findfirst(x->x==8, FreeDOF) #vector location of load

#define pedestrian load at node 5
#the vertical DOF number for node 5 is 14
#Apply load starting from 10 sec to 80 sec
tFirst = 0:0.01:9.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=10:0.01:80 #seconds    Turn on the load
ptSecond=a1*sin.(ω1*tSecond)
tThird=80.01:0.01:130         #Turn off the load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt2Spline = Spline1D(tAll,ptAll)
pt2DOF=findfirst(x->x==14, FreeDOF)  #vector location of load

#define pedestrian load at node 7
#the vertical DOF number for node 7 is 20
#Apply load starting from 20 sec to 90 sec
tFirst = 0:0.01:19.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=20:0.01:90 #seconds    Turn on the load
ptSecond=a1*sin.(ω1*tSecond)
tThird=90.01:0.01:130          #Turn off the load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt3Spline = Spline1D(tAll,ptAll)
pt3DOF=findfirst(x->x==20, FreeDOF)  #vector location of load

#define pedestrian load at node 9
#the vertical DOF number for node 9 is 26
#Apply load starting from 30 sec to 100 sec
tFirst = 0:0.01:29.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=30:0.01:100 #seconds    Turn on load
ptSecond=a1*sin.(ω1*tSecond)
tThird=100.01:0.01:130          #Turn off load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt4Spline = Spline1D(tAll,ptAll)
pt4DOF=findfirst(x->x==26, FreeDOF)  #vector location of load

#define pedestrian load at node 11
#the vertical DOF number for node 11 is 32
#Apply load starting from 40 sec to 110 sec
tFirst = 0:0.01:39.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=40:0.01:110 #seconds    Turn on load
ptSecond=a1*sin.(ω1*tSecond)
tThird=110.01:0.01:130          #Turn off load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt5Spline = Spline1D(tAll,ptAll)
pt5DOF=findfirst(x->x==32, FreeDOF)  #vector location of load


#define pedestrian load at node 13
#the vertical DOF number for node 13 is 39
#Apply load starting from 50 sec to 120 sec
tFirst = 0:0.01:49.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=50:0.01:120 #seconds    Turn on load
ptSecond=a1*sin.(ω1*tSecond)
tThird=120.01:0.01:130          #Turn off load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt6Spline = Spline1D(tAll,ptAll)
pt6DOF=findfirst(x->x==38, FreeDOF) #vector location of load

#define pedestrian load at node 15
#the vertical DOF number for node 15 is 44
#Apply load starting from 60 sec to 130 sec
tFirst = 0:0.01:59.9 #seconds
ptFirst=zeros(length(tFirst))
tSecond=60:0.01:130 #seconds    Turn on load
ptSecond=a1*sin.(ω1*tSecond)
tThird=130.01:0.01:130.5        #Turn off load
ptThird=zeros(length(tThird))
tAll=[tFirst;tSecond;tThird]
ptAll=[ptFirst;ptSecond;ptThird]
pt7Spline = Spline1D(tAll,ptAll)
pt7DOF=findfirst(x->x==44, FreeDOF) #vector location of load
# using Plots
# plot(tAll,pt2Spline(tAll))
# plot(tAll,ptAll)

#
# pt=zeros(size(M)[1])
#     pt[pt2DOF]=pt2Spline(60.0)


#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt1DOF, pt2Spline, pt2DOF, pt3Spline, pt3DOF, pt4Spline, pt4DOF, pt5Spline, pt5DOF, pt6Spline, pt6DOF, pt7Spline, pt7DOF = p

    pt=zeros(size(M)[1])  #define vector of all the pts at each free DOF

    pt[pt1DOF]=pt1Spline(t[1])  #assign load to node 3, vertical DOF
    pt[pt2DOF]=pt2Spline(t[1])  #assign load to node 5, vertical DOF
    pt[pt3DOF]=pt3Spline(t[1])  #assign load to node 7, vertical DOF
    pt[pt4DOF]=pt4Spline(t[1])  #assign load to node 9, vertical DOF
    pt[pt5DOF]=pt5Spline(t[1])  #assign load to node 11, vertical DOF
    pt[pt6DOF]=pt6Spline(t[1])  #assign load to node 13, vertical DOF
    pt[pt7DOF]=pt7Spline(t[1])  #assign load to node 15, vertical DOF
    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1]

end
 #                                    u_dot0       u_0         trange
  prob = SecondOrderODEProblem(mdof, zeros(size(M)[1]), zeros(size(M)[1]), (0.,200.),(M, C, K, pt1Spline, pt1DOF, pt2Spline, pt2DOF, pt3Spline, pt3DOF, pt4Spline, pt4DOF, pt5Spline, pt5DOF, pt6Spline, pt6DOF, pt7Spline, pt7DOF))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:200)

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
plot(t,u[:,8],xlabel="t[sec]",ylabel="displacement[in]",label="3")
plot!(t,u[:,14],xlabel="t[sec]",ylabel="displacement[in]",label="5")
plot!(t,u[:,20],xlabel="t[sec]",ylabel="displacement[in]",label="7")
plot!(t,u[:,26],xlabel="t[sec]",ylabel="displacement[in]",label="9")  #plot displacement at DOF 26
plot!(t,u[:,32],xlabel="t[sec]",ylabel="displacement[in]",label="11")
plot!(t,u[:,38],xlabel="t[sec]",ylabel="displacement[in]",label="13")
plot!(t,u[:,44],xlabel="t[sec]",ylabel="displacement[in]",label="15")
