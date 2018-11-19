function [D, r_ab] = GPSGlobal(xa_pred,P_pred,Phi,R, hg, a)
%This function calculates the Global position of the agent using a GPS measurement

i = 2*a+1;
j = 2*a+2;

Pa = P_pred(i:j,i:j);
Phia = Phi(i:j,i:j);


h_pred = sqrt((xa_pred(1))^2+(xa_pred(2))^2);

%calculate residual range

ra = hg-h_pred;

% define jacobians of the measurement function step K-1
rxa = -(xa_pred(1))/h_pred;
rya = -(xa_pred(2))/h_pred;

Ha = [rxa rya];

% Hab = [Ha,Hb];
Sa = R + Ha*Pa*Ha';
temp = sqrt(1/Sa);

Da0 = -(Phia^(-1)*Pa*Ha')*temp;

r_ab = temp*ra;

D = zeros(6,1);

D(i:j) = Da0;

end

