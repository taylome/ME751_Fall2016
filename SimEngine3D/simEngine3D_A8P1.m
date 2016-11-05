%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Oct 2016
%
% Homework 8 - Problem #1
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

clc
close all

[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\me751_HW07.mdl'];

Model_IVD = simEngine3D(model_name);

[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\me751_HW08P1.mdl'];

Model_Dyn = simEngine3D(model_name);

%%
close all

figure();
subplot(2,1,1);
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(1,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(2,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Force (N)');
title(['ME751 - Reaction Forces from the Revolute Joint on the Pendulum in the Global Frame - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('FX','FY','FZ');
set(gca(),'FontSize',16)
subplot(2,1,2);
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(1,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(2,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Force (N)');
title('ME751 - Reaction Forces from the Revolute Joint on the Pendulum in the Global Frame - Inverse Dynamics');
legend('FX','FY','FZ');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1);
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(4,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(5,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.ConstraintReactions(1).Body2(6,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Torque (Nm)');
title(['ME751 - Reaction Torques from the Revolute Joint on the Pendulum in the Global Frame - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('TX','TY','TZ');
set(gca(),'FontSize',16)
ylim([-250,250]);
subplot(2,1,2);
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(4,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(5,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.ConstraintReactions(1).Body2(6,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Torque (Nm)');
title('ME751 - Reaction Torques from the Revolute Joint on the Pendulum in the Global Frame - Inverse Dynamics');
legend('TX','TY','TZ');
set(gca(),'FontSize',16)
ylim([-250,250]);

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(1,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(2,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Displacement (m)');
title(['ME751 - Pendulum CG Postion in the Global Frame - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).q(1,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).q(2,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).q(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Displacement (m)');
title('ME751 - Pendulum CG Postion in the Global Frame - Inverse Dynamics');
legend('X','Y','Z');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(4,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(5,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(6,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).q(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters');
title(['ME751 - Pendulum Euler Parameters - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).q(4,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).q(5,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).q(6,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).q(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters');
title('ME751 - Pendulum Euler Parameters - Inverse Dynamics');
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(1,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(2,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Velocity (m)');
title(['ME751 - Pendulum CG Velocity in the Global Frame - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).qd(1,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).qd(2,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Velocity (m)');
title('ME751 - Pendulum CG Velocity in the Global Frame - Inverse Dynamics');
legend('X','Y','Z');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(4,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(5,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(6,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qd(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters Dot');
title(['ME751 - Pendulum Euler Parameters Dot - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).qd(4,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).qd(5,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qd(6,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qd(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters');
title('ME751 - Pendulum Euler Parameters Dot - Inverse Dynamics');
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(1,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(2,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Acceleration (m/s)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('X','Y','Z');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(1,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(2,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Acceleration (m/s)');
title('ME751 - Pendulum CG Accel in the Global Frame - Inverse Dynamics');
legend('X','Y','Z');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(4,:),'linewidth',3);
hold on
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(5,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(6,:),'linewidth',3);
plot(Model_Dyn.time,Model_Dyn.bodies(1).qdd(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters Dot');
title(['ME751 - Pendulum Euler Parameters DDot - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(4,:),'linewidth',3);
hold on
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(5,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(6,:),'linewidth',3);
plot(Model_IVD.time,Model_IVD.bodies(1).qdd(7,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Euler Parameters');
title('ME751 - Pendulum Euler Parameters DDot - Inverse Dynamics');
legend('e0','e1','e2','e3');
set(gca(),'FontSize',16)

figure();
subplot(2,1,1)
plot(Model_Dyn.time,Model_Dyn.NR);
xlabel('Time (s)');
ylabel('Number of NR Iterations');
title(['ME751 - Number of NR Iterations - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
set(gca(),'FontSize',16)
subplot(2,1,2)
plot(Model_Dyn.time,Model_Dyn.CondPsi);
xlabel('Time (s)');
ylabel('Number of NR Iterations');
title(['ME751 - Final Condition Number of Psi - Dynamics - Time Step: ',num2str(Model_Dyn.simulation.stepSize),'s']);
set(gca(),'FontSize',16)
