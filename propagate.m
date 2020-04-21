function P_dot = propagate(P,A,C,R,Q)
%PROPAGATE Summary of this function goes here
%   Detailed explanation goes here
P = reshape(P,size(A));
P_dot = A*P +P*A'+Q-P*C'*R^-1*C*P;
P_dot = P_dot(:);
end

