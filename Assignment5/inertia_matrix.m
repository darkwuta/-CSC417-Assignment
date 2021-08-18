syms phi1 phi2 x_0 x_1 x_2 y_0 y_1 y_2 z_0 z_1 z_2;
phi2max = @(phi1)1-phi1;
xxx = @(phi1,phi2)((1-phi1-phi2)*x_0+phi1*x_1+phi2*x_2)^3;
yyy = @(phi1,phi2)((1-phi1-phi2)*y_0+phi1*y_1+phi2*y_2)^3;
zzz = @(phi1,phi2)((1-phi1-phi2)*z_0+phi1*z_1+phi2*z_2)^3;
xxy = @(phi1,phi2)((1-phi1-phi2)*x_0+phi1*x_1+phi2*x_2)^2*((1-phi1-phi2)*y_0+phi1*y_1+phi2*y_2);
yyz = @(phi1,phi2)((1-phi1-phi2)*y_0+phi1*y_1+phi2*y_2)^2*((1-phi1-phi2)*z_0+phi1*z_1+phi2*z_2);
xxz = @(phi1,phi2)((1-phi1-phi2)*x_0+phi1*x_1+phi2*x_2)^2*((1-phi1-phi2)*z_0+phi1*z_1+phi2*z_2);

q1 = int(int(xxx,phi2,0,phi2max),phi1,0,1);
q2 = int(int(yyy,phi2,0,phi2max),phi1,0,1);
q3 = int(int(zzz,phi2,0,phi2max),phi1,0,1);

q4 = int(int(xxy,phi2,0,phi2max),phi1,0,1);
q5 = int(int(yyz,phi2,0,phi2max),phi1,0,1);
q6 = int(int(xxz,phi2,0,phi2max),phi1,0,1);

q1
q2
q3
q4
q5
q6

