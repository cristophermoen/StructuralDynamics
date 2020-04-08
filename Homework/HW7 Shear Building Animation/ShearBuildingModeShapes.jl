using LinearAlgebra   #for eig function
using Makie   #to show mode shapes

#show vibration mode shapes for a three story shear building

#floors
m1 = 8000/2.281  #kg   #first floor, 20 psf dead load
m2 = 4000/2.281  #kg   #second floor (or roof), 10 psf dead load
m3 = 2000/2.281  #kg   #third floor (or roof), 5 psf dead load

#shear walls per story
E=2E+11/8  #N/m^2   #assume concrete
b=3     #m  width
t=0.2   #m  thickness
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


#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
wn_squared=eigvals(K,M)

#solve for eigenvectors K*phi_n=wn^2*M*phi_n
modeshapes=eigvecs(K,M)

L=h

#show modes

#pick mode to display
DisplayMode=1

#define the mode
u1=modeshapes[1,DisplayMode]
u2=modeshapes[2,DisplayMode]
u3=modeshapes[3,DisplayMode]
mode=[u1, u2, u3]

#normalize the mode
VectorScaleIndex=argmax(abs.([u1, u2, u3]))
mode=mode./mode[VectorScaleIndex]

#scale the mode
ModeDisplayScale=1.0
mode=mode.*ModeDisplayScale


#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum(mode)-4.0
Ymin=0-4.0
Xmax=maximum(mode)+4.0
Ymax=3*L+4.0
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)


#define initial position
XYPos1=Node(mode[1])
XYPos2=Node(mode[2])
XYPos3=Node(mode[3])


#draw walls

#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L], [-L/2+xy, L], [-L/2, 0], [L/2, 0]], XYPos1), color = :lightgray, show_axis = false, limits=limits, overdraw =:true)
#second floor
poly!(scene1, lift((xy1,xy2)->Point2f0[[L/2+xy2, 2*L], [-L/2+xy2, 2*L], [-L/2+xy1, L], [L/2+xy1, L]], XYPos1,XYPos2), color = :lightgray, show_axis = false, limits=limits)
#third floor
poly!(scene1, lift((xy2,xy3)->Point2f0[[L/2+xy3, 3*L], [-L/2+xy3, 3*L], [-L/2+xy2, 2*L], [L/2+xy2, 2*L]], XYPos2, XYPos3), color = :lightgray, show_axis = false, limits=limits)



#draw floors
#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L+L/16], [-L/2+xy, L+L/16], [-L/2+xy, L-L/16], [L/2+xy, L-L/16]], XYPos1), color = :darkgray, show_axis = false, limits=limits)
#second floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, 2*L+L/16], [-L/2+xy, 2*L+L/16], [-L/2+xy, 2*L-L/16], [L/2+xy, 2*L-L/16]], XYPos2), color = :darkgray, show_axis = false, limits=limits)
#third floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, 3*L+L/16], [-L/2+xy, 3*L+L/16], [-L/2+xy, 3*L-L/16], [L/2+xy, 3*L-L/16]], XYPos3), color = :darkgray, show_axis = false, limits=limits)


#draw ground and building centerline
lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)

#draw masses

scatter!(scene1, lift(xy->Point2f0[(xy,L)],XYPos1), marker = [:circle], limits=limits, color= :red, markersize=1)  #first floor mass
scatter!(scene1, lift(xy->Point2f0[(xy,2*L)],XYPos2), marker = [:circle], limits=limits, color= :red, markersize=0.5)  #second floor mass
scatter!(scene1, lift(xy->Point2f0[(xy,3*L)],XYPos3), marker = [:circle], limits=limits, color= :red, markersize=0.25)  #third floor mass

#add mode number, natural frequency, and natural period
wn=sqrt(wn_squared[DisplayMode])
fn=round(wn/(2*pi),digits=1)
Tn=round(1/fn, digits=3)

title(scene1,"Mode $DisplayMode, fn=$fn cycles/sec., Tn=$Tn sec.")
