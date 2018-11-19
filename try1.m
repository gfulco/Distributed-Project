clc
clear all
close all

%% agent 1 model

dt = 0.1;
R = 0.002; % say measurement variance is the same for all the agents (same sensors)
ag1out = zeros(10000,3);

x10 = [2 2 2]';
A1 = eye(3);
B1 = eye(3)*dt;
u1 = [1 1 1]';
P10 = ones(3,3)/100;
Q1 = eye(3,3)/200;
%P1_12 = zeros(3,3);
%P1_21 = P1_12';
Phi10 = eye(3);

%% first propagation for both agent

%agent 1
P1_pred = A1*P10*A1'+Q1;
x1_pred = A1*x10+B1*u1;
Phi1 = A1*Phi10;

%% first update (no message)

%agent 1
x1 = x1_pred;
P1 = P1_pred;
%P1_12 = P1_21;

x11 = A1*x10+B1*u1;
ag1out(1,:) = x11'-x1';

for k=2:10000
%% propagation step K
clc
%agent 1
P1_pred = A1*P1*A1'+Q1
x1_pred = A1*x1+B1*u1
Phi1 = A1*Phi1;

%% suppose agent1 takes an absolute measurement of its position (range to the origin)

x11 = A1*x11+B1*u1;

h10 = sqrt((x11(1))^2+(x11(2))^2+(x11(3))^2)+(2*rand()-1)*0.002;

%moreover use the measurement function with the propagation stage datas

h10_pred = sqrt((x1_pred(1))^2+(x1_pred(2))^2+(x1_pred(3))^2);

%calculate residual range

r10 = h10-h10_pred

% define jacobians of the measurement function step K-1

rx10 = (x1_pred(1))/h10_pred;
ry10 = (x1_pred(2))/h10_pred;
rz10 = (x1_pred(3))/h10_pred;

H10 = [rx10 ry10 rz10];

S10 = R + H10*P1_pred*H10'

D10 = (Phi1^(-1)*P1_pred*H10')*S10^(-1/2)

r_10 = S10^(-1/2)*r10

%% update stage with agent 1 that measured its absolute position

% agent 1 update

x1 = x1_pred+Phi1*D10*r_10;
P1 = P1_pred-Phi1*(D10*D10')*Phi1';

ag1out(k,:) = x11'-x1';
end

figure(1)
clf
plot(ag1out);





