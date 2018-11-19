function [xhat,P] = UpdateStep( xhat_pred,P_pred,Phi,D,r,a)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
i = 2*a+1;
j = 2*a+2;

xhat = xhat_pred+Phi(i:j,i:j)*D(i:j)*r;
P = zeros(6,6);


for k = 1:2:size(xhat_pred)-2
    for l = k:2:size(xhat_pred)
       if k == l && k ~= i
           continue
       elseif k == l && k == i
           P(i:j,i:j) = P_pred-Phi(i:j,i:j)*(D(i:j)*D(i:j)')*Phi(i:j,i:j)';
       elseif k ~= l 
            P(k:k+1,l:l+1) = P_pred(k:k+1,l:l+1) - D(k:k+1)*D(l,l+1);
            P(l:l+1,k:k+1) = (P_pred(k:k+1,l:l+1) - D(k:k+1)*D(l,l+1))';
       end
    end


   
end