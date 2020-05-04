using LinearAlgebra  #for calculating eigenvalues
using Plots
using ArbNumerics
import SpecialFunctions: besselj, bessely


T=55    #lb
ζ=0.05  #damping ratio
ro=50   #lb/m
c=sqrt(T/ro)
R=1
u=zeros(21,1000)
rt=[2.4 5.5 8.65 11.79]

for i=1:21
	for j=1:1000
		for k=1:4
			a=besselj(2,rt[k])
			b=besselj(1,rt[k])^2*rt[k]^2
			c=4/R^(2)*a/b
			global u[i,j]=u[i,j]+c*cos(2*rt[k]*j/10)*besselj(0,rt[k]*sqrt((i-11)^2)/R/10)*exp(-ζ*c*rt[k]*j/(R*10))
		end
	end
end

#p1=plot(1:10, u[:,1])
#p2=plot(1:10,u[:,23])

@gif for i=1:200
	plot(1:21,u[:,i],ylim=[-1,1])
end every 1
