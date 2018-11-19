function [xhat_pred,P_pred,Phi_next] = PropagationStep(xhat,A,B,u,P,Q,Phi,index)
%Propagation Step of the DCL Algorithm, always the same for every agent,
%with measurement or not.

i = 2*index+1;
j = 2*index+2;
P_pred = P;
P_pred(i:j,i:j) = A*P(i:j,i:j)*A'+Q;
xhat_pred = A*xhat(i:j)+B*u;
Phi_next = Phi;
Phi_next(i:j,i:j) = A*Phi(i:j,i:j);

end

