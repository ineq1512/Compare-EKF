function P_dot = hybridpropagate(P,A,L,Q)
%PROPAGATE Summary of this function goes here
%   Detailed explanation goes here
P = reshape(P,size(A));
P_dot = A*P +P*A'+L*Q*L';
P_dot = P_dot(:);
end