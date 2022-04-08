function [ x ] = robot( x0,u,h )
%% State space system of 2-wheel mobile robot
%% physical parameters and design parameters 
mc = 30; % mass of robot - [kg]
mw = 1; % mass of a wheel - [kg]
r = 0.15; % radius of wheels - [m]
b = 0.75; % length between wheels - [m]
d = 0.3; % length between wheel center and gravity center - [m]
Ic = 15.625; % moment of inertia of mobile robot - [kg*m^2]
Iw = 0.005;%  the wheel with a motor about the wheel axis - [kg*m^2]
Im = 0.0025; % the wheel with a motor about the wheel diameter - [kg*m^2]
%%
I = mc*d^2 + 2*mw*b^2 + Ic + 2*Im;
m = mc + 2*mw;
th = x0(3);
S = [ 0.5*r*cos(th), 0.5*r*cos(th);
	  0.5*r*sin(th), 0.5*r*sin(th);
	  0.5*r/b, -0.5*r/b];
M = [ r^2/(4*b^2) * (m*b^2 + I) + Iw, r^2/(4*b^2) * (m*b^2 - I);
      r^2/(4*b^2) * (m*b^2 - I),      r^2/(4*b^2) * (m*b^2 + I) + Iw];
B = [ 1, 0, 0.5*cos(th), 0.5*sin(th);
	  0, 1, 0.5*cos(th), 0.5*sin(th)];
    
eta_dot = h* S*x0(4:5);
dth = eta_dot(3);
V = [ 0, r^2/(2*b) * mc*d*dth;
      -r^2/(2*b) * mc*d*dth, 0];

nu_dot = h* (M)\( B*u - V*x0(4:5));
nu_ = x0(4:5) + nu_dot;
eta_ = x0(1:3) + eta_dot;
x = [eta_;nu_];
end

