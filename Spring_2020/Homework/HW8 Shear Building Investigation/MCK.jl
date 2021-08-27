#Calculate the dynamic response of a three story shear building

#floors
m1 = 2400/2.281  #kg   #first floor
m2 = 3200/2.281  #kg   #second floor (or roof)
m3 = 1000/2.281  #kg   #third floor (or roof)

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
