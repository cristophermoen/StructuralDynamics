using LinearAlgebra   #for eig function
using Makie   #to show mode shapes

#Calculate the free vibration response of a two story shear building

#floors
m1 = 8000/2.281  #kg   #first floor, 20 psf dead load
m2 = 4000/2.281  #kg   #second floor (or roof), 10 psf dead load

#columns
E=2E+11  #N/m^2
Ic=723/12^4/3.281^4  #m^4   #W14x68 steel column, strong axis
L=13/3.281  #m  column height
k1 = 12*E*4*Ic/L^3
k2 = 12*E*2*Ic/L^3

#define mass, stiffness matrices
M = [m1 0
     0 m2]
K = [k1+k2 -k2
     -k2   k2]


#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
wn_squared=eigvals(K,M)

#solve for eigenvectors K*phi_n=wn^2*M*phi_n
modeshapes=eigvecs(K,M)


#show modes

#pick mode to display
DisplayMode=2

#define the mode
u1=modeshapes[1,DisplayMode]
u2=modeshapes[2,DisplayMode]
mode=[u1, u2]

#normalize the mode
VectorScaleIndex=argmax(abs.([u1, u2]))
mode=mode./mode[VectorScaleIndex]

#scale the mode
ModeDisplayScale=1.0
mode=mode.*ModeDisplayScale


#define window and limits
scene1=Scene(resolution=(1000,1000))
Xmin=minimum(mode)-4.0
Ymin=0-4.0
Xmax=maximum(mode)+4.0
Ymax=2*L+4.0
Xrange=Xmax-Xmin
Yrange=Ymax-Ymin
limits = FRect(Xmin, Ymin, Xrange, Yrange)


#define initial position
XYPos1=Node(mode[1])
XYPos2=Node(mode[2])


#draw floors
#first floor
poly!(scene1, lift(xy->Point2f0[[L/2+xy, L+L/16], [-L/2+xy, L+L/16], [-L/2+xy, L-L/16], [L/2+xy, L-L/16]], XYPos1), color = :lightgray, show_axis = false, limits=limits)
#second floor
poly!(scene1, lift((xy1,xy2)->Point2f0[[L/2+xy1+xy2, 2*L+L/16], [-L/2+xy1+xy2, 2*L+L/16], [-L/2+xy1+xy2, 2*L-L/16], [L/2+xy1+xy2, 2*L-L/16]], XYPos1, XYPos2), color = :lightgray, show_axis = false, limits=limits)

#draw ground and building centerline
lines!([0, 0], [Ymin, Ymax], color= :black, linestyle= :dashdot)
lines!([Xmin, Xmax], [0, 0], color= :green, lineweight=6)

#draw masses

scatter!(scene1, lift(xy->Point2f0[(xy,L)],XYPos1), marker = [:circle], limits=limits, color= :red, markersize=1)  #first floor mass
scatter!(scene1, lift((xy1,xy2)->Point2f0[(xy1+xy2,2*L)],XYPos1, XYPos2), marker = [:circle], limits=limits, color= :red, markersize=0.5)  #second floor mass


#draw columns
#use shape function and coordinate system from
#https://www.youtube.com/watch?v=K7lfyAldj5k

function BeamShape(q1,q2,q3,q4,L,x, offset)

    a0=q1
    a1=q2
    a2=1/L^2*(-3*q1-2*q2*L+3*q3-q4*L)
    a3=1/L^3*(2*q1+q2*L-2*q3+q4*L)

    w=a0 .+a1.*x .+a2.*x.^2 .+a3.*x.^3 .+offset

end


x=0:L/10:L
q1=0
q2=0
q4=0

lines!(lift(xy->BeamShape(q1,q2,xy,q4,L,x, L/2),XYPos1),x, linewidth=12)   #first level right column
lines!(lift(xy->BeamShape(q1,q2,xy,q4,L,x, -L/2),XYPos1),x, linewidth=12) #first level left column

lines!(lift((xy1, xy2)->BeamShape(q1,q2,xy2,q4,L,x,L/2+xy1),XYPos1,XYPos2),x .+L, linewidth=6) #second level right column
lines!(lift((xy1, xy2)->BeamShape(q1,q2,xy2,q4,L,x,-L/2+xy1),XYPos1,XYPos2),x .+L, linewidth=6) #second level left column

#add mode number, natural frequency, and natural period
wn=sqrt(wn_squared[DisplayMode])
fn=round(wn/(2*pi),digits=1)
Tn=round(1/fn, digits=3)

title(scene1,"Mode $DisplayMode, fn=$fn cycles/sec., Tn=$Tn sec.")
