using OrdinaryDiffEq  #ODE solver
using Dierckx  #for interpolating an aribtrary loading
using LinearAlgebra  #for calculating eigenvalues
using Plots
#floors
m1 = 1000000  #kg   #first floor
m2 = 1000000
m3 = 1000000
m4 = 0.02*(m1+m2+m3)
#shear walls per story
E=2E+11/8  #N/m^2   #assume concrete
b=3     #m  width
t=0.1   #m  thickness
A=b*t   #top view cross-sectional area of shear wall
ν=0.20
G=E/(2*(1+ν))
h=13/3.281  #m  story height
k1 = 2*G*A/h
k2 = 2*G*A/h
k3 = G*A/h
#determine natural frequency
M0 = [m1 0 0
     0 m2 0
     0 0  m3]

K0 = [k1+k2 -k2 0
     -k2   k2+k3 -k3
     0     -k3    k3]
ωn_squared=eigvals(K0,M0)
ωn=sqrt.(ωn_squared)
k4=m4*ωn_squared[1]
#define mass, stiffness, and damping matrices
M = [m1 0 0   0
     0 m2 0   0
     0 0  m3  0
     0 0  0  m4]

K = [k1+k2 -k2   0   0
     -k2   k2+k3 -k3 0
     0     -k3   k3+k4  -k4
     0     0     -k4   k4]
w=0.00:0.01:2.1*ωn[1]
u1max=zeros(length(w))
u0max=zeros(length(w))
a=zeros(length(w))
b=zeros(length(w))
for i = 1:length(w)
    a[i]=w[i]/ωn[1]
    t1 = 0:0.01:30 #seconds
    a1 = 100000
    ω1 = w[i]
    pt1 = a1*sin.(ω1*t1)
    pt1Spline = Spline1D(t1,pt1)

        #at DOF 2
        t2 = 0:0.01:30 #seconds
        a2 = 0
        ω2 = 320
        pt2 = a2*cos.(ω2*t2)
        pt2Spline = Spline1D(t2,pt2)
        #at DOF 3
        t3 = 0:0.01:30 #seconds
        a3 = 0
        ω3 = 320
        pt3 = a3*cos.(ω3*t3)
        pt3Spline = Spline1D(t3,pt3)

        #at DOF 4
        t4 = 0:0.01:30 #seconds
        a4 = 0.0
        ω4 = 1.0
        pt4 = a4*cos.(ω4*t4)
        pt4Spline = Spline1D(t4,pt4)
        function damper(ddu, du, u, p, t)

            M, K, pt1Spline, pt2Spline, pt3Spline, pt4Spline = p
            pt1=pt1Spline(t[1])
            pt2=pt2Spline(t[1])
            pt3=pt3Spline(t[1])
            pt4=pt4Spline(t[1])
            pt=[pt1; pt2; pt3; pt4]

            ddu[:,1] =- inv(M)*K*u[:,1] + inv(M)*pt[:,1]

        end
         #                                    u_dot0            u_0     trange
          prob = SecondOrderODEProblem(damper, [0.; 0.; 0.; 0.], [0.0; 0.0; 0.0; 0.0], (0.,10.),(M, K, pt1Spline, pt2Spline, pt3Spline, pt4Spline))
          sol = solve(prob, DPRKN6(),tstops=0:0.01:10)
          u1_dot=(x->x[1]).(sol.u)
          u2_dot=(x->x[2]).(sol.u)
          u3_dot=(x->x[3]).(sol.u)
          u4_dot=(x->x[4]).(sol.u)
          u1=(x->x[5]).(sol.u)
          u2=(x->x[6]).(sol.u)
          u3=(x->x[7]).(sol.u)
          u4=(x->x[8]).(sol.u)
          t=sol.t
          u1max[i]=maximum(u3)

          function mdof(ddu, du, u, p, t)

              M0, K0, pt1Spline, pt2Spline, pt3Spline = p

              pt1=pt1Spline(t[1])
              pt2=pt2Spline(t[1])
              pt3=pt3Spline(t[1])
              pt=[pt1; pt2; pt3]

              ddu[:,1] = - inv(M0)*K0*u[:,1] + inv(M0)*pt[:,1]

          end
           #                                    u_dot0            u_0     trange
            prob = SecondOrderODEProblem(mdof, [0.; 0.; 0.], [0.0; 0.0; 0.0], (0.,10.),(M0, K0, pt1Spline, pt2Spline, pt3Spline))
            sol = solve(prob, DPRKN6(),tstops=0:0.01:10)

            u1_dot=(x->x[1]).(sol.u)
            u2_dot=(x->x[2]).(sol.u)
            u3_dot=(x->x[3]).(sol.u)
            u01=(x->x[4]).(sol.u)
            u02=(x->x[5]).(sol.u)
            u03=(x->x[6]).(sol.u)
            t=sol.t
            u0max[i]=maximum(u03)
            b[i]=u1max[i]./u0max[i]

  end
plot(a,b,xlabel="forcing frequency/natural frequency",ylabel="u1/u0",title="3 stories building")
