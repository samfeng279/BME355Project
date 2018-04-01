% Simply call 
%
%       >> double_pendulum_init
%
% to run the double pendulum simulation with the below parameters. This
% script calls double_pendulum.
%
%   ---------------------------------------------------------------------

phi1                = 7*pi/4;
dtphi1              = pi/12;
phi2                = 2*pi+0.6;
dtphi2              = 0;
g                   = 9.81; 
m1                  = 3.285; 
m2                  = 1.095; 
l1                  = 0.4318; 
l2                  = 0.07136;
duration            = .5;
fps                 = 100;
movie               = true;

clc; figure;

restLengthTA = StabilityModel.tibialisLength(pi/2);
TA = HillTypeMuscle(2000, .6*restLengthTA, .4*restLengthTA);
restLengthS = StabilityModel.soleusLength(pi/2);
S = HillTypeMuscle(16000, .6*restLengthS, .4*restLengthS);

interval=[0, duration];
ivp=[phi1; dtphi1; phi2; dtphi2; g; m1; m2; l1; l2;0;0];

double_pendulum(ivp, duration, fps, movie, TA, S);