function [ D_up ] = UpdateMessage( P, D,Kb,Ka,a,b )
%UpdateMessage sent by interim master to all the agents in the network

i = 2*a+1;
j = 2*b+1;

D_up = zeros(6,1);
for m = 1:2:size(P)
    if m == 2*a+1 || m == 2*b+1
        D_up(m:m+1) = D(m:m+1);
    else
        D_up(m:m+1) = P(m:m+1,j:j+1)*Kb-P(m:m+1,i:i+1)*Ka;
    end


end

