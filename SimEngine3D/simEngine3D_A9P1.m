%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Nov 2016
%
% Homework 9 - Problem #1
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

clc
close all

[pathstr,~,~] = fileparts(mfilename('fullpath'));

model_name = [pathstr,'\models\me751_HW09P1_FullNewton.mdl'];
Model_FullNewton = simEngine3D(model_name);

model_name = [pathstr,'\models\me751_HW09P1_ModifiedNewton.mdl'];
Model_ModifiedNewton = simEngine3D(model_name);

model_name = [pathstr,'\models\me751_HW09P1_QuasiModifiedNewton.mdl'];
Model_QuasiModifiedNewton = simEngine3D(model_name);


figure();
hold on
plot(Model_FullNewton.time,Model_FullNewton.NR);
plot(Model_ModifiedNewton.time,Model_ModifiedNewton.NR);
plot(Model_QuasiModifiedNewton.time,Model_QuasiModifiedNewton.NR);
xlabel('Time (s)');
ylabel('Number of NR Iterations');
title(['ME751 - Number of NR Iterations - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s']);
legend('Full Newton','Modified Newton','Quasi Modified Newton');
set(gca(),'FontSize',16)

figure();
subplot(3,1,1);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(1,:)-Model_QuasiModifiedNewton.bodies(1).qdd(1,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Quasi Modified NR"']);
legend('X');
set(gca(),'FontSize',16)
subplot(3,1,2);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(2,:)-Model_QuasiModifiedNewton.bodies(1).qdd(2,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Quasi Modified NR"']);
legend('Y');
set(gca(),'FontSize',16)
subplot(3,1,3);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(3,:)-Model_QuasiModifiedNewton.bodies(1).qdd(3,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Quasi Modified NR"']);
legend('Z');
set(gca(),'FontSize',16)

figure();
subplot(3,1,1);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(1,:)-Model_ModifiedNewton.bodies(1).qdd(1,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Modified NR"']);
legend('X');
set(gca(),'FontSize',16)
subplot(3,1,2);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(2,:)-Model_ModifiedNewton.bodies(1).qdd(2,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Modified NR"']);
legend('Y');
set(gca(),'FontSize',16)
subplot(3,1,3);
plot(Model_FullNewton.time,Model_FullNewton.bodies(1).qdd(3,:)-Model_ModifiedNewton.bodies(1).qdd(3,:));
grid on
xlabel('Time (s)');
ylabel('Accel Difference (m/s^2)');
title(['ME751 - Pendulum CG Accel in the Global Frame - Dynamics - Time Step: ',num2str(Model_FullNewton.simulation.stepSize),'s "Full NR" - "Modified NR"']);
legend('Z');
set(gca(),'FontSize',16)
