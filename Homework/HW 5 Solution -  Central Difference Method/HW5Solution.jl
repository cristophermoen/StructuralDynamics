#function style before code is debugged
# *******************

#sample problem
#time
dt=0.01
totalt=5.0
t=0.0:dt:totalt   #define time range

#define loading function
pt=zeros(length(t))

#dynamic system properties
m=300
k=100000

ζ=0.05

#initial conditions
uo=0.1
duo=0.0



ωn=sqrt(k/m)
fn=ωn/(2*pi)

#critical damping ratio
ccr=2*m*ωn
c=ζ*ccr

#define acceleration at t=0
dduo=-c/m*duo-k/m*uo+pt[1]/m

#define ghost node displacement at t=0
uMinus1=uo-duo*dt+dduo*dt^2/2

#calculate kHat
kHat=m/dt^2+c/(2*dt)

a=m/dt^2-c/(2*dt)
b=k-(2*m)/dt^2

for i=1:(length(t)-1)

    #first step you need uMinus1

    pHat[i]=pt[i]-a*u[i-1]-b*u[i]

    u[i+1]=pHat[i]/kHat

end

# *******************






#function style after code is debugged
# *******************
function CentralDifference(dt,totalt,pt,m,k,ζ,uo,duo)

    ωn=sqrt(k/m)
    fn=ωn/(2*pi)

    #critical damping ratio
    ccr=2*m*ωn
    c=ζ*ccr

    #define acceleration at t=0
    dduo=-c/m*duo-k/m*uo+pt[1]/m

    #define ghost node displacement at t=0
    uMinus1=uo-duo*dt+dduo*dt^2/2

    #calculate kHat
    kHat=m/dt^2+c/(2*dt)

    a=m/dt^2-c/(2*dt)
    b=k-(2*m)/dt^2

    for i=1:(length(t)-1)

        #first step you need uMinus1

        pHat[i]=pt[i]-a*u[i-1]-b*u[i]

        u[i+1]=pHat[i]/kHat

    end

    return u

end


#sample problem
#time
dt=0.01
totalt=5.0
t=0.0:dt:totalt   #define time range

#define loading function
pt=zeros(length(t))

#dynamic system properties
m=300
k=100000

ζ=0.05

#initial conditions
uo=0.1
duo=0.0

u=CentralDifference(dt,totalt,pt,m,k,ζ,uo,duo)


# *******************
