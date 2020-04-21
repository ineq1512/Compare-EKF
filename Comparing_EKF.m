% Program build to compare perfomance between Continuous - time EKF,
% hybrid EKF and the discrete - time EKF.
% The dynamic model is reffered to Example 13.2 in "Dan Simon - Optimal
% state estimation_ Kalman, H-infinity, and nonlinear approaches-Wiley-Interscience (2006)"

%%
clear
close all
clc
%% Parameter
p0 = 0.0034; % lb-sec^2/ft^4
g  = 32.2;   % ft/sec^2
k  = 22000;   % ft
%% Sampling time and simulation time
Ts = 0.4e-3; Tsim = 16; Ns = Tsim/Ts;
%% Parameter and initial states
x0    = [100000 -6000 1/2000]';
xhat0 = [100000 -6100 1/2500]';
P0    = [500 0       0;
         0 20000     0;
         0 0     1/250000;];

varw = [0 0 0]';
varv = 100;
mode = 3; % 1 = continuous EKF, 2 = hybrid EKF, 3 = discrete EKF 

% Tm is the sampling time for filter, it equal to Ts in continuous time
% case and equal to 0.5s in other cases
if mode ==1
    Tm =Ts; 
else
    Tm = 0.5;
end
%% State variable 
x = zeros(3,Ns+1);
x(:,1) = x0;
xhat = zeros(3,Tsim/Tm+1);
xhat(:,1) = xhat0;
y = zeros(1,Ns+1);
%% Filter' parameter
P2 = zeros(3,3);
P1 = P0;  
C  = [1 0 0];
H  = C;
L  = eye(3);
Q  = L*diag(varw)*L';
M  = 1;
R  = varv;

%% Simulation
tic
for i =1:1:(Tsim/Ts)
    
    %Dynamic model
        w = [wgn(1,1,varw(1),'linear'); wgn(1,1,varw(2),'linear'); wgn(1,1,varw(3),'linear')];
        v = wgn(1,1,varv,'linear');
        [~,x_ode] =ode45(@(t,x_ode) dynamicmodel(x_ode,w), [(i-1)*Ts, i*Ts],x(:,i));
        x(:,i+1) = x_ode(size(x_ode,1),:)';
        y(i+1) = x(1,i+1)+v; % update the measurement 
        
    % State estimation
    if mod(i,Tm/Ts) == 0
        if mode ==1
            iest = i;
        else 
            iest = fix(i*Ts/Tm);
        end
        %Update the state matrix A
        A21 = (-p0*exp(-xhat(1,iest)/k)*xhat(2,iest)^2)/(2*k*xhat(3,iest));
        A22 = (p0*exp(-xhat(1,iest)/k)*xhat(2,iest))/xhat(3,iest);
        A23 = (-p0*exp(-xhat(1,iest)/k)*xhat(2,iest)^2)/(2*xhat(3,iest)^2);
        A = [0   1   0;
             A21 A22 A23;
             0    0  0;];
         % Update the estimate state xhat
        [xhat(:,iest+1),P2] = estimate(xhat(:,iest),y(i+1),A,P1,w,varw,varv,[0, Tm],mode); 
  
        P1 = P2; % Update the matrix P

    end

end
toc
%% Plot
d = fix(2*Tm/Ts)
figure(1)
idx =0:Tm:Tsim;
plot(idx,xhat(2,:));
figure(2)
if mode ==1
    plot(idx,xhat(2,:)-x(2,:));
else
    plot(idx,xhat(2,:)-x(2,idx*d+1));
end
% hold on 
% plot(k,y);