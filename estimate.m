function [state_out,P2] = estimate(state_in,y,A,P1,w,varw,varv,tspan, mode)
%ESTIMATE Summary of this function goes here
%   Detailed explanation goes here

L = eye(3);
Q = diag(varw);
R = varv;
C = [1 0 0];
H = C;
F=c2d(A,zeros(3,1),0.5);
if mode == 1 % Continuous EKF
    K = P1*C'*R^-1;
    [~,state_ode] = ode45(@(t,state_ode) dynamicest(state_ode,y,w,K), tspan,state_in); % propagate xhat
    state_out = state_ode(size(state_ode,1),:)';
    [~,P_ode] = ode45(@(t,P_ode) propagate(P_ode,A,C,R,Q), tspan, P1); % propagate P
    P2 = P_ode(size(P_ode,1),:)';
    P2 = reshape(P2,size(A));
elseif mode == 2 % hybrid EKF
    % Integrate the state estimate and its covariance from time k-1 to k
    [~,state_ode] = ode45(@(t,state_ode) dynamicmodel(state_ode,w), tspan,state_in); % propagate xhat
    state_out = state_ode(size(state_ode,1),:)';
    [~,P_ode] = ode45(@(t,P_ode) hybridpropagate(P_ode,A,L,Q), tspan, P1); % propagate P-
    P2 = P_ode(size(P_ode,1),:)';
    P2 = reshape(P2,size(A));
    K = P2*H'*((H*P2*H'+R)^-1);
    state_out = state_out + (K*(y-state_out(1)));
    P2 = (eye(1)-K*H)*P2*(eye(1)-K*H)'+K*R*K';
elseif mode == 3 % discrete EKF
    P2 = F*P1*F' +L*Q*L';
    state_out = F*state_in;
    K = P2*H'*(H*P2*H'+R)^-1;
    state_out = state_out +K*(y-state_out(1));
    P2 = (eye(1)- K*H)*P2*(eye(1)- K*H)'+K*R*K';
end
end

