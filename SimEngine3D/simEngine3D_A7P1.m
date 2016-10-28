%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Oct 2016
%
% Homework 7 - Problem #1
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

clc
close all

[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\me751_HW07.mdl'];

Model = simEngine3D(model_name);

close all

figure()
plot(Model.time,Model.ConstraintReactions(1).Body2(1,:),'linewidth',3);
hold on
plot(Model.time,Model.ConstraintReactions(1).Body2(2,:),'linewidth',3);
plot(Model.time,Model.ConstraintReactions(1).Body2(3,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Force (N)');
title('ME751 - Reaction Forces from the Revolute Joint on the Pendulum in the Global Frame');
legend('FX','FY','FZ');
set(gca(),'FontSize',16)


figure()
plot(Model.time,Model.ConstraintReactions(1).Body2(4,:),'linewidth',3);
hold on
plot(Model.time,Model.ConstraintReactions(1).Body2(5,:),'linewidth',3);
plot(Model.time,Model.ConstraintReactions(1).Body2(6,:),'linewidth',3);
grid on
xlabel('Time (s)');
ylabel('Torque (Nm)');
title('ME751 - Reaction Torques from the Revolute Joint on the Pendulum in the Global Frame');
legend('TX','TY','TZ');
set(gca(),'FontSize',16)
