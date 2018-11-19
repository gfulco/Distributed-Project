function [ Pab, Phiab ] = LandmarkMessage(Pa,Pb,Phia, Phib, a, b )
% The interim Master gets a message with xb_hat, Pb_pred, Phib and updates
% its matrices. No need to save xb_hat cause it's already available to
% each agent in this implementation.
 
i = 2*a+1;
j = 2*a+2;
k = 2*b+1;
l = 2*b+2;

Pab = Pa;
Phiab = Phia;

Pab(k:l,k:l) = Pb(k:l,k:l);
Phiab(k:l,k:l) = Phib(k:l,k:l);


end

