%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Nov 2016
%
% Homework 8 - Problem #2
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

% %Pendulum 1 Initial Position& Mass Props
% L = 2
% m = 78
% 1/12*m*(0.1^2+0.1^2)
% 1/12*m*(0.1^2+(2*L)^2)
%
% qa.e0 = cos((pi()/2)/2);
% e = [0;1;0]*sin((pi()/2)/2);
% qa.e1 = e(1);
% qa.e2 = e(2);
% qa.e3 = e(3);
% 
% qb.e0 = cos((pi()/2)/2);
% e = [0;0;1]*sin((pi()/2)/2);
% qb.e1 = e(1);
% qb.e2 = e(2);
% qb.e3 = e(3);
% 
% q = zeros(4,1);
% q(1) = qa.e0 * qb.e0 - qa.e1 * qb.e1 - qa.e2 * qb.e2 - qa.e3 * qb.e3;
% q(2) = qa.e0 * qb.e1 + qa.e1 * qb.e0 - qa.e3 * qb.e2 + qa.e2 * qb.e3;
% q(3) = qa.e0 * qb.e2 + qa.e2 * qb.e0 + qa.e3 * qb.e1 - qa.e1 * qb.e3;
% q(4) = qa.e0 * qb.e3 + qa.e3 * qb.e0 - qa.e2 * qb.e1 + qa.e1 * qb.e2;
% 
% r = [0;2*sin(pi()/2);-2*cos(pi()/2)];
% q = [r;...
%      q]
 

% %Pendulum 2 Initial Position & Mass Props
% L = 2/2
% m = 78/2
% 1/12*m*(0.1^2+0.1^2)
% 1/12*m*(0.1^2+(2*L)^2)
% 
% qa.e0 = cos((pi()/2)/2);
% e = [0;1;0]*sin((pi()/2)/2);
% qa.e1 = e(1);
% qa.e2 = e(2);
% qa.e3 = e(3);
% 
% qb.e0 = cos((pi())/2);
% e = [0;0;1]*sin((pi())/2);
% qb.e1 = e(1);
% qb.e2 = e(2);
% qb.e3 = e(3);
% 
% q = zeros(4,1);
% q(1) = qa.e0 * qb.e0 - qa.e1 * qb.e1 - qa.e2 * qb.e2 - qa.e3 * qb.e3;
% q(2) = qa.e0 * qb.e1 + qa.e1 * qb.e0 - qa.e3 * qb.e2 + qa.e2 * qb.e3;
% q(3) = qa.e0 * qb.e2 + qa.e2 * qb.e0 + qa.e3 * qb.e1 - qa.e1 * qb.e3;
% q(4) = qa.e0 * qb.e3 + qa.e3 * qb.e0 - qa.e2 * qb.e1 + qa.e1 * qb.e2;
% 
% r = [0;2;-1];
% q = [r;...
%      q]
 
 
clc
close all


[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\me751_HW08P2.mdl'];

Model = simEngine3D(model_name);

figure();
plot(Model.time,Model.bodies(1).q(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(1).q(2,:),'linewidth',3);
plot(Model.time,Model.bodies(1).q(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Displacement (m)');
title(['ME751 - Pendulum 1 CG Postion in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)


figure();
plot(Model.time,Model.bodies(1).qd(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(1).qd(2,:),'linewidth',3);
plot(Model.time,Model.bodies(1).qd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Velocity (m)');
title(['ME751 - Pendulum 1 CG Velocity in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)



figure();
plot(Model.time,Model.bodies(1).qdd(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(1).qdd(2,:),'linewidth',3);
plot(Model.time,Model.bodies(1).qdd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Acceleration (m/s)');
title(['ME751 - Pendulum 1 CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)



figure();
plot(Model.time,Model.bodies(2).q(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(2).q(2,:),'linewidth',3);
plot(Model.time,Model.bodies(2).q(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Displacement (m)');
title(['ME751 - Pendulum 2 CG Postion in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)


figure();
plot(Model.time,Model.bodies(2).qd(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(2).qd(2,:),'linewidth',3);
plot(Model.time,Model.bodies(2).qd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Velocity (m)');
title(['ME751 - Pendulum 2 CG Velocity in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)


figure();
plot(Model.time,Model.bodies(2).qdd(1,:),'linewidth',3);
hold on
plot(Model.time,Model.bodies(2).qdd(2,:),'linewidth',3);
plot(Model.time,Model.bodies(2).qdd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Acceleration (m/s)');
title(['ME751 - Pendulum 2 CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)


figure();
plot(Model.bodies(1).q(2,:),Model.bodies(1).q(3,:),'linewidth',3);
hold on
plot(Model.bodies(2).q(2,:),Model.bodies(2).q(3,:),'linewidth',3);
grid on
xlabel('Pendulum CG Y Position');
ylabel('Pendulum CG Z Position');
title(['ME751 - Pendulum CG Position in the Global Frame - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
legend('Pendulum 1','Pendulum 2');
set(gca(),'FontSize',16)
axis image


Phi2_norm = zeros(size(Model.Phi,2),1);
for i = 1:length(Phi2_norm)
    Phi2_norm(i) = norm(Model.Phi(6:10,i));
end


figure();
subplot(2,1,1);
plot(Model.time,Phi2_norm,'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Norm2 of Revolute 2 Constraints');
title(['ME751 - Norm2 of the Second Revolute Position Constraint Equations - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
set(gca(),'FontSize',16)
subplot(2,1,2);
plot(Model.time,Model.NR,'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Iteration Count');
title(['ME751 - NR Iteration Count - Dynamics - Time Step: ',num2str(Model.simulation.stepSize),'s']);
set(gca(),'FontSize',16)

