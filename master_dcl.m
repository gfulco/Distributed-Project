%% restart

clc
clear all
close all

%% agent 1 model

dt = .1;
R = 0.002; % say measurement variance is the same for all the agents (same sensors)


xhat10 = [0 0 0]';
A1 = eye(3);
B1 = eye(3)*dt;
u1 = [1 1 1]';
P10 = ones(3,3)/200;
Q1 = eye(3,3)/100;
P1_12 = zeros(3,3);
P1_21 = P1_12';
Phi10 = eye(3);

%% agent 2 model

xhat20 = [2 2 2]';
A2 = eye(3);
B2 = eye(3)*dt;
u2 = [1 1 1]';
P20 = ones(3,3)/200;
Q2 = eye(3,3)/100;
P2_12 = zeros(3,3);
P2_21 = P2_12';
Phi20 = eye(3);

%% first propagation for both agent

%agent 1
P1_pred = A1*P10*A1'+Q1;
xhat1_pred = A1*xhat10+B1*u1;
Phi1 = A1*Phi10;

%agent 2
P2_pred = A2*P20*A2'+Q2;
xhat2_pred = A2*xhat20+B2*u2;
Phi2 = A2*Phi20;


%% first update (no message)

%agent 1
xhat1 = xhat1_pred;
P1 = P1_pred;
P1_12 = P1_12;
P1_21 = P1_12';


%agent 2
xhat2 = xhat2_pred;
P2 = P2_pred;
P2_12 = P2_12;
P2_21 = P2_12';


xhat11 = A1*xhat10+B1*u1;
xhat22 = A2*xhat20+B2*u2;


ag1out = zeros(2000,3);
ag2out = zeros(2000,3);


for(k=2:20000)

%% k-th iteration propagation 

%agent 1
P1_pred = A1*P1*A1'+Q1;
xhat1_pred = A1*xhat1+B1*u1;
Phi1 = A1*Phi1;

%agent 2
P2_pred = A2*P2*A2'+Q2;
xhat2_pred = A2*xhat2+B2*u2;
Phi2 = A2*Phi2;


%% suppose a message is received in the k-th iteration
% now agent 1 know datas from propagation stage of agent 2

%say the measurement taken by agent 1 is the real distance from agent 2
%plus an error +/- 0.2

h = sqrt((xhat1(1)-xhat2(1))^2+(xhat1(2)-xhat2(2))^2+(xhat1(3)-xhat2(3))^2)+(2*rand()-1)*0.002;

%moreover use the measurement function with the propagation stage datas

h_pred = sqrt((xhat1_pred(1)-xhat2_pred(1))^2+(xhat1_pred(2)-xhat2_pred(2))^2+(xhat1_pred(3)-xhat2_pred(3))^2);

%calculate residual range

r12 = h-h_pred;

% define jacobians of the measurement function step K-1
rxhat1 = -(xhat1_pred(1)-xhat2_pred(1))/h_pred;
rxhat2 = -(xhat1_pred(1)-xhat2_pred(1))/h_pred;
ry1 = -(xhat1_pred(2)-xhat2_pred(2))/h_pred;
ry2 = -(xhat1_pred(2)-xhat2_pred(2))/h_pred;
rz1 = -(xhat1_pred(3)-xhat2_pred(3))/h_pred;
rz2 = -(xhat1_pred(3)-xhat2_pred(3))/h_pred;

H1 = [rxhat1 ry1 rz1];
H2 = [rxhat2 ry2 rz2];

H12 = [H1,H2];

S12 = R + H1*P1_pred*H1' + H2*P2_pred*H2'-H1*Phi1*P1_12*Phi2'*H2' - H2*Phi2*P1_21*Phi1'*H1';

D1 = (Phi1^(-1)*Phi1*P1_12*Phi2'*H2'-Phi1^(-1)*P1_pred*H1')*S12^(-1/2);
D2 = (Phi2^(-1)*P2_pred*H2'-P1_21*Phi1'*H1')*S12^(-1/2);

r_12 = S12^(-1/2)*r12;

%% update stage with agent 2 that receive the response

% agent 1 update

xhat1 = xhat1_pred+Phi1*D1*r_12;
P1 = P1_pred-Phi1*(D1*D1')*Phi1';
P1_12 = P1_12-D1*D2';
P1_21 = P1_12';

clc
[xhat1_pred Phi1*D1*r_12]
pause(1e-3)

%agent 2 update

xhat2 = xhat2_pred+Phi2*D2*r_12;
P2 = P2_pred-Phi2*(D2*D2')*Phi2';
P2_12 = P2_12-D1*D2';
P2_21 = P2_12';

xhat11 = A1*xhat11+B1*u1;
xhat22 = A2*xhat22+B2*u2;

ag1out(k,:) = xhat11'-xhat1';
ag2out(k,:) = xhat22'-xhat2';
end

figure(1)
clf
subplot(211)
plot(ag1out);
subplot(212)
plot(ag2out);
