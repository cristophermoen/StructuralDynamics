using LinearAlgebra   #for eig function
# using Makie   #to show mode shapes

#show vibration mode shapes for a three story shear building

#floors
m1 = 1200000  #kg   #first floor, 20 psf dead load
m2 = 1200000  #kg   #second floor (or roof), 10 psf dead load
m3 = 600000  #kg   #third floor (or roof), 5 psf dead load

#shear walls per story
k1 = 63045606000
k2 = 63045606000
k3 = 63045606000

#define mass, stiffness, and damping matrices
M = [m1 0 0
     0 m2 0
     0 0  m3]

K = [k1+k2 -k2 0
     -k2   k2+k3 -k3
     0     -k3    k3]


#solve for eigenvalues of K*phi_n=wn^2*M*phi_n
wn_squared=eigvals(K,M)

wn=sqrt.(wn_squared)

#solve for eigenvectors K*phi_n=wn^2*M*phi_n
modeshapes=eigvecs(K,M)

L=15/3.281

#show modes

#pick mode to display
DisplayMode=1

#define the mode
u1=modeshapes[1,DisplayMode]
u2=modeshapes[2,DisplayMode]
u3=modeshapes[3,DisplayMode]
mode=[u1, u2, u3]

u0=[0;0;0.5]

ϕn1=modeshapes[:,1]
VectorScaleIndex=argmax(abs.([ϕn1[1], ϕn1[2], ϕn1[3]]))
ϕn1=ϕn1./ϕn1[VectorScaleIndex]

Mn1=transpose(ϕn1)*M*ϕn1

qn10=transpose(ϕn1)*M*u0/Mn1

t=0.0:0.01:3.0


qn1=qn10.*cos.(wn[1]*t)


Mode1Contribution=zeros(3,length(t))

for i=1:length(t)
    Mode1Contribution[:,i]=qn1[i]*ϕn1
end

# Mode2Contribution=zeros(3,length(t))
# for i=1:length(t)
#     Mode2Contribution[:,i]=qn2[i]*ϕn2
# end
#
# Mode3Contribution=zeros(3,length(t))
# for i=1:length(t)
#     Mode3Contribution[:,i]=qn3[i]*ϕn3
# end
#
# #compare this to HW6
# Total=Mode1Contribution+Mode2Contribution+Mode3Contribution

using Plots
#certain floor

plot(t,Mode1Contribution[1,:])


# plot!(t,Mode2Contribution)
# plot!(t,Mode3Contribution)


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
