function [D, r_ab,Kb,Ka] = KalmanGain(xa_pred,xb_pred,P_pred,Phi,R, h, a, b )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

i = 2*a+1;
j = 2*a+2;
k = 2*b+1;
l = 2*b+2;

Pa = P_pred(i:j,i:j);
Pb = P_pred(k:l,k:l);
Pab = P_pred(i:j,k:l);
Phia = Phi(i:j,i:j);
Phib = Phi(k:l,k:l);

h_pred = sqrt((xb_pred(1)-xa_pred(1))^2+(xb_pred(2)-xa_pred(2))^2);

%calculate residual range

rab = h-h_pred;

% define jacobians of the measurement function step K-1
rxa = (xb_pred(1)-xa_pred(1))/h_pred;
rxb = (xb_pred(1)-xa_pred(1))/h_pred;
rya = (xb_pred(2)-xa_pred(2))/h_pred;
ryb = (xb_pred(2)-xa_pred(2))/h_pred;

Ha = [rxa rya];
Hb = [rxb ryb];

% Hab = [Ha,Hb];

Sab = R + Ha*Pa*Ha' + Hb*Pb*Hb' -Ha*Phia*Pab*Phib'*Hb' - Hb*Phib*Pab'*Phia'*Ha';

Da = (Phia^(-1)*Phia*Pab*Phib'*Hb'-Phia^(-1)*Pa*Ha').*sqrt(1/Sab);
Db = (Phib^(-1)*Pb*Hb'-Pab'*Phia'*Ha')*Sab^(-1/2);

r_ab = Sab^(-1/2)*rab;


D = zeros(6,1);

D(i:j) = Da;
D(k:l) = Db;

Kb = Phi(k:l,k:l)'*Hb'*Sab^(-1/2);
Ka = Phi(i:j,i:j)'*Ha'*Sab^(-1/2);

end

