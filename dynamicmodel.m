function [state_dot] = dynamicmodel(state,noise)
%PMSM Summary of this function goes here
% Dynamic model

%% Parameter
p0 = 0.0034; % lb-sec^2/ft^4
g  = 32.2;   % ft/sec^2
k  = 22000;   % ft
%% System
state_dot = zeros(3,1);
x1 = state(1);
x2 = state(2);
x3 = state(3);
w1 = noise(1);
w2 = noise(2);
w3 = noise(3);

%% Dynamic model
state_dot(1) = x2+w1;
state_dot(2) = ((p0*exp(-x1/k)*x2^2)/(2*x3))-g+w2;
state_dot(3) = w3;
end

