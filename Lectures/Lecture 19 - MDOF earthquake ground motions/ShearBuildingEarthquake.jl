using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading
using LinearAlgebra  #for calculating eigenvalues
using DelimitedFiles #use for importing earthquake acceleration file



cd(@__DIR__)  #change working directory to that of this file

#read in ground acceleration data, note dims=(rows,columns) in text file that has 711 lines and 8 columns
#Imperial Valley 1979, Imperial County Center Grounds ground motion
data = readdlm("IV1979acc.txt",Float64,dims=(711,8))

#rearrange ground accleration data as a vector
data=transpose(data)   #transpose matrix
ddu_g=reshape(data,(5688,1))  #reshape matrix into a vector, note 8*711=5688
ddu_g=vec(ddu_g)   #tell Julia to make the data a 1D vector

#units come in as cm/sec^2 from strongmotion.org, let's change them to m/sec^2
ddu_g=ddu_g/100

#define time range of the earthquake assuming the first acc reading occurs at t=0 seconds
t_eq=collect(range(0,length=5688,stop=5687*0.01))   #total time is 5687*0.01 seconds

#define ground acceleration as a spline for later interpolation
EQSpline = Spline1D(t_eq,ddu_g)


#Calculate the dynamic response of a three story shear building

#floors
m1 = 240000/2.281  #kg
m2 = 320000/2.281  #kg
m3 = 100000/2.281  #kg

#shear walls per story
E=2E+11/8  #N/m^2   #assume concrete
b=3     #m  width
t=0.1   #m  thickness
A=b*t   #top view cross-sectional area of shear wall
ν=0.20
G=E/(2*(1+ν))
h=13/3.281  #m  story height
k1 = G*A/h
k2 = G*A/h
k3 = G*A/h

#define mass, stiffness, and damping matrices
M = [m1 0 0
     0 m2 0
     0 0  m3]

K = [k1+k2 -k2 0
     -k2   k2+k3 -k3
     0     -k3    k3]

#calculate natural frequencies
#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
ωn_squared=eigvals(K,M)
ωn=sqrt.(ωn_squared)

#Raleigh damping
ωi=ωn[1]   #set first frequency anchor point i
ωj=ωn[3]   #set second frequency anchor point j
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


#earthquake influence vector ι
#define which Free DOF feel the ground acceleration
#in this example it is all three horizontal DOF
#if there were vertical DOF, they would not feel the horizontal ground acceleration
ι=[1;1;1]

#write equation of motion function
function mdof(ddu, du, u, p, t)

    M, C, K, pt1Spline, pt2Spline, pt3Spline, EQSpline, ι= p

    pt1=pt1Spline(t[1])
    pt2=pt2Spline(t[1])
    pt3=pt3Spline(t[1])
    pt=[pt1; pt2; pt3]

    ddu_g=EQSpline(t[1])

    ddu[:,1] = -inv(M)*C*du[:,1] - inv(M)*K*u[:,1] + inv(M)*pt[:,1] - ι*ddu_g

end
 #                                    u_dot0      u_0     trange
  prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.], [0.0; 0.0; 0.0], (0.,50.),(M, C, K, pt1Spline, pt2Spline, pt3Spline, EQSpline, ι))
  sol = solve(prob, DPRKN6(),tstops=0:0.01:50)

u1_dot=(x->x[1]).(sol.u)
u2_dot=(x->x[2]).(sol.u)
u3_dot=(x->x[3]).(sol.u)
u1=(x->x[4]).(sol.u)
u2=(x->x[5]).(sol.u)
u3=(x->x[6]).(sol.u)

t=sol.t
dt=t_eq[2]-t_eq[1]

#approximate ground displacement integrating ddu_g twice
du_g=cumsum(diff(ddu_g)*dt) #ground velocity
u_g=[0; 0; cumsum(diff(du_g)*dt)]  #ground displacement


function ShearBuildingAnimation(u1, u2, u3, u_g, L, DeformationScale)

    #animate the three story shear building

    #Scale deformation to study response history
    u1=u1*DefScale
    u2=u2*DefScale
    u3=u3*DefScale
    u_g=u_g*DefScale

    #define window and limits
    scene1=Scene(resolution=(1000,1000))
    Xmin=minimum([u1;u2;u3])./DeformationScale-4.0
    Ymin=0-4.0
    Xmax=maximum([u1;u2;u3])./DeformationScale+4.0
    Ymax=3*L+4.0
    Xrange=Xmax-Xmin
    Yrange=Ymax-Ymin
    limits = FRect(Xmin, Ymin, Xrange, Yrange)

    #define initial position
    XYPos1=Node(u1[1])
    XYPos2=Node(u2[1])
    XYPos3=Node(u3[1])
    XYPos4=Node(u_g[1])

    #draw walls
    #first floor
    poly!(scene1, lift((xy, xy4)->Point2f0[[L/2+xy+xy4, L], [-L/2+xy+xy4, L], [-L/2+xy4, 0], [L/2+xy4, 0]], XYPos1, XYPos4), color = :lightgray, show_axis = false, limits=limits, overdraw =:true)
    #second floor
    poly!(scene1, lift((xy1,xy2,xy4)->Point2f0[[L/2+xy2+xy4, 2*L], [-L/2+xy2+xy4, 2*L], [-L/2+xy1+xy4, L], [L/2+xy1+xy4, L]], XYPos1,XYPos2, XYPos4), color = :lightgray, show_axis = false, limits=limits)
    #third floor
    poly!(scene1, lift((xy2,xy3,xy4)->Point2f0[[L/2+xy3+xy4, 3*L], [-L/2+xy3+xy4, 3*L], [-L/2+xy2+xy4, 2*L], [L/2+xy2+xy4, 2*L]], XYPos2, XYPos3, XYPos4), color = :lightgray, show_axis = false, limits=limits)

    #draw floors
    #first floor
    poly!(scene1, lift((xy,xy4)->Point2f0[[L/2+xy+xy4, L+L/16], [-L/2+xy+xy4, L+L/16], [-L/2+xy+xy4, L-L/16], [L/2+xy+xy4, L-L/16]], XYPos1, XYPos4), color = :darkgray, show_axis = false, limits=limits)
    #second floor
    poly!(scene1, lift((xy,xy4)->Point2f0[[L/2+xy+xy4, 2*L+L/16], [-L/2+xy+xy4, 2*L+L/16], [-L/2+xy+xy4, 2*L-L/16], [L/2+xy+xy4, 2*L-L/16]], XYPos2, XYPos4), color = :darkgray, show_axis = false, limits=limits)
    #third floor
    poly!(scene1, lift((xy,xy4)->Point2f0[[L/2+xy+xy4, 3*L+L/16], [-L/2+xy+xy4, 3*L+L/16], [-L/2+xy+xy4, 3*L-L/16], [L/2+xy+xy4, 3*L-L/16]], XYPos3, XYPos4), color = :darkgray, show_axis = false, limits=limits)

    #draw building centerline
    lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)

    #draw ground
    poly!(scene1, lift(xy4->Point2f0[[3*L/4+xy4, 0], [-3*L/4+xy4, 0], [-3*L/4+xy4, -L/2], [3*L/4+xy4, -L/2]], XYPos4), color = :green, show_axis = false, limits=limits)

    #draw masses
    scatter!(scene1, lift((xy,xy4)->Point2f0[(xy+xy4,L)],XYPos1, XYPos4), marker = [:circle], limits=limits, color= :red, markersize=m1/100000)  #first floor mass
    scatter!(scene1, lift((xy, xy4)->Point2f0[(xy+xy4,2*L)],XYPos2, XYPos4), marker = [:circle], limits=limits, color= :red, markersize=m2/100000)  #second floor mass
    scatter!(scene1, lift((xy, xy4)->Point2f0[(xy+xy4,3*L)],XYPos3, XYPos4), marker = [:circle], limits=limits, color= :red, markersize=m3/100000)  #third floor mass

    #animate
    record(scene1, "animation.mp4", range(1, stop = length(u1), step=5)) do i

        if i==1
            sleep(1)
        end

        XYPos1[]=(u1[i])
        XYPos2[]=(u2[i])
        XYPos3[]=(u3[i])
        XYPos4[]=(u_g[i])

        # sleep(1/24)

    end

end


# **********
using Makie
DeformationScale=1000.0
ShearBuildingAnimation(u1, u2, u3, u_g, L, DeformationScale)
