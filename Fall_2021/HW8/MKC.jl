#Define the structural properties of the mystery building.

#floors
m1 = 40000  #kg   #first floor
m2 = 20000  #kg   #second floor 
m3 = 10000  #kg   #third floor 

k1 =  5000000  #N/m
k2 = 10000000  #N/m
k3 = 10000000  #N/m

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

#' Calculate the modal natural frequencies.
ωn = sqrt.(ωn_squared)

#' Calculate the mode shape natural frequencies in Hz.
fn = ωn ./(2 * π)

#' Calculate the mode shape natural periods in seconds.
Tn = 1 ./ fn

#Raleigh damping
ωi=ωn[1]   #set first frequency anchor point i
ωj=ωn[3]   #set second frequency anchor point j
ζi=0.02    #modal viscous damping ratio in mode i
ζj=0.10    #modal viscous damping ratio in mode j

a0,a1=2*inv([1/ωi ωi;1/ωj ωj])*[ζi;ζj]  #Eq. 11.4.9 from Chopra
C=a0*M+a1*K
