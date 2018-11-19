%% restart

clc
clear all
close all


%% datas

dt = 0.1;
iter = 20000;
R = 0.5; % say measurement variance is the same for all the agents (same sensors)
Rg = 0.05;
maxrange = 3;
x0 = [0, 0, 1, 1, -1, -1]';
x = zeros(6,iter);
x(:,1) = x0; % real agents position
A = eye(2);
B = eye(2)*dt;
D0 = zeros(6,1);
agentsout = zeros(6,iter); 
u = [1 1]';

%% each agents keep the system model for each agent also they keep the covariance matrix and the Phi matrix for each agent int the network
%agent 1 model

xhat = zeros(6,iter);
xhat(:,1) = [50,50,100,100,-50,-50]';

P1 = zeros(6,6);
P1(1:2,1:2) = [5, 0; 0,5];
Q1 = eye(2)/20;
Phi1 = zeros(6,6);
Phi1(1:2,1:2) = eye(2,2); 
D1 = zeros(6,1);
xhat1 = zeros(2,iter);
xhat1(:,1) = zeros(2,1); % predicted agent position
% xhat1_pred = zeros(2,1);
% P1_pred = zeros(6,6);

% agent 2 model


P2 = zeros(6,6);
P2(3:4,3:4) = [5, 0; 0, 5];
Q2 = eye(2)/20;
Phi2 = zeros(6,6);
Phi2(3:4,3:4) = eye(2,2); 
D2 = zeros(6,1);
xhat2 = zeros(2,iter);
xhat2(:,1) = zeros(2,1); % predicted agent position
% xhat2_pred = zeros(2,1);
% P2_pred = zeros(6,6);

% agent 3 model

P3 = zeros(6,6);
P3(5:6,5:6) = [5, 0; 0, 5];
Q3 = eye(2)/20;
Phi3 = zeros(6,6);
Phi3(5:6,5:6) = eye(2,2); 
D3 = zeros(6,1);
xhat3 = zeros(2,iter);
xhat3(:,1) = zeros(2,1); % predicted agent position
% xhat3_pred = zeros(2,1);
% P3_pred = zeros(6,6);

% Create a matrix that contain the error of the filter at each step.

agentsout(:,1) = x(:,1) - xhat(:,1);

%% Say the algorithm is launched and no measurement is taken in the first step k = 1

k=1;


u1 = u.*(2*rand(2,1)-1)*2;
u2 = u.*(2*rand(2,1)-1)*5;
u3 = u.*(2*rand(2,1)-1)*3;

%at every step the "real" position of the agent is updated.
x(:,k+1) = eye(6)*x0+eye(6)*dt*[u1;u2;u3];

%Propagation step
[xhat1_pred,P1_pred,Phi1] = PropagationStep(xhat(:,k),A,B,u1,P1,Q1,Phi1,0);
[xhat2_pred,P2_pred,Phi2] = PropagationStep(xhat(:,k),A,B,u2,P2,Q2,Phi2,1);
[xhat3_pred,P3_pred,Phi3] = PropagationStep(xhat(:,k),A,B,u3,P3,Q3,Phi3,2);

%no measurement, we pass only the predicted matrices to the function,
%others are all zero so we pass D0 as default D and r = 0 and THE
%SAME INDEX a b

[xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D0,0,0);
[xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D0,0,1);
[xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D0,0,2);

% error updating 

xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);

%% successive iteration k = 2:iter-1



for k = 2:iter
    
    u1 = u.*(2*rand(2,1)-1)*2;
    u2 = u.*(2*rand(2,1)-1)*5;
    u3 = u.*(2*rand(2,1)-1)*3;
    
    %at every step the "real" position of the agent is updated.
    x(:,k+1) = eye(6)*x(:,k)+eye(6)*dt*[u1;u2;u3];
    
    %Propagation step
    [xhat1_pred,P1_pred,Phi1] = PropagationStep(xhat(:,k),A,B,u1,P1,Q1,Phi1,0);
    [xhat2_pred,P2_pred,Phi2] = PropagationStep(xhat(:,k),A,B,u2,P2,Q2,Phi2,1);
    [xhat3_pred,P3_pred,Phi3] = PropagationStep(xhat(:,k),A,B,u3,P3,Q3,Phi3,2);
    
    % check if there measurement are taken
    
    h(1) = MeasureDistance(x(1:2,k+1),x(3:4,k+1),maxrange);
    h(2) = MeasureDistance(x(1:2,k+1),x(5:6,k+1),maxrange);
    h(3) = MeasureDistance(x(3:4,k+1),x(5:6,k+1),maxrange);
    
    

    if mod(k,7) == 0
        h1 = MeasureDistance(x(1:2,k+1),[0 0],NaN);
        [D1,r1] = GPSGlobal(xhat1_pred,P1_pred,Phi1,Rg, h1, 0);
        [xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D1,r1,0);
        [xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D1,r1,1);
        [xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D1,r1,2);
        xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
        agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);
    continue
    elseif ~isnan(h(1)) 
             [P1_pred,Phi1] = LandmarkMessage(P1_pred,P2_pred,Phi1,Phi2,0,1);
             [D1,r_12,K2,K1] = KalmanGain(xhat1_pred,xhat2_pred,P1_pred,Phi1,R, h(1), 0, 1);
             D1 = UpdateMessage(P1_pred,D1,K2,K1,0,1);
             D2 = UpdateMessage(P2_pred,D1,K2,K1,0,1);
             D3 = UpdateMessage(P3_pred,D1,K2,K1,0,1);
             [xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D1,r_12,0);
             [xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D2,r_12,1);
             [xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D3,r_12,2);
             xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
             agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);
             continue
    elseif ~isnan(h(2))
            [P1_pred,Phi1] = LandmarkMessage(P1_pred,P3_pred,Phi1,Phi3,0,2);
            [D1,r_13,K3,K1] = KalmanGain(xhat1_pred,xhat3_pred,P1_pred,Phi1,R, h(2), 0, 2);
            D1 = UpdateMessage(P1_pred,D1,K3,K1,0,2);
            D2 = UpdateMessage(P2_pred,D1,K3,K1,0,2);
            D3 = UpdateMessage(P3_pred,D1,K3,K1,0,2);
            [xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D1,r_13,0);
            [xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D2,r_13,1);
            [xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D3,r_13,2);
            xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
            agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);
            continue
    elseif ~isnan(h(3)) 
            [P2_pred,Phi2] = LandmarkMessage(P2_pred,P3_pred,Phi2,Phi3,1,2);
            [D2,r_23,K3,K2] = KalmanGain(xhat2_pred,xhat3_pred,P2_pred,Phi2,R, h(3), 1, 2);
            D1 = UpdateMessage(P1_pred,D2,K3,K2,1,2);
            D2 = UpdateMessage(P2_pred,D2,K3,K2,1,2);
            D3 = UpdateMessage(P3_pred,D2,K3,K2,1,2);
            [xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D1,r_23,0);
            [xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D2,r_23,1);
            [xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D3,r_23,2);
            xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
            agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);
            continue
   else
            [xhat1(:,k+1),P1] = UpdateStep(xhat1_pred,P1_pred,Phi1,D0,0,0); 
            [xhat2(:,k+1),P2] = UpdateStep(xhat2_pred,P2_pred,Phi2,D0,0,1);
            [xhat3(:,k+1),P3] = UpdateStep(xhat3_pred,P3_pred,Phi3,D0,0,2);
            xhat(:,k+1) = [xhat1(:,k+1);xhat2(:,k+1);xhat3(:,k+1)];
            agentsout(:,k+1) = x(:,k+1) - xhat(:,k+1);
 
     end
end


%% plotting 

figure('Name','Agent 1')
subplot(311)
plot(x(1,:),x(2,:))
title('Real Position')
subplot(312)
plot(xhat(1,:),xhat(2,:))
title('Estimated Position')
subplot(313)
plot(1:size(agentsout,2),agentsout(1:2,:)')
title('Estimation Error')

figure('Name','Agent 2')
subplot(311)
plot(x(3,:),x(4,:))
title('Real Position')
subplot(312)
plot(xhat(3,:),xhat(4,:))
title('Estimated Position')
subplot(313)
plot(1:size(agentsout,2),agentsout(3:4,:)')
title('Estimation Error')

figure('Name','Agent 3')
subplot(311)
plot(x(5,:),x(6,:))
title('Real Position')
subplot(312)
plot(xhat(5,:),xhat(6,:))
title('Estimated Position')
subplot(313)
plot(1:size(agentsout,2),agentsout(5:6,:)')
title('Estimation Error')
