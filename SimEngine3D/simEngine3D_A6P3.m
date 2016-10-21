%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Oct 2016
%
% Homework 6 - Problem #3
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

function Model = simEngine3D_A6P3()
clc
close all

[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\me751_HW06.mdl'];

Model = simEngine3D(model_name);

%generate Plots
for i = 1:size(Model.bodies(1).points,1)
    pos = zeros(3,size(Model.bodies(1).q,2));
    vel = zeros(3,size(Model.bodies(1).qd,2));
    acc = zeros(3,size(Model.bodies(1).qdd,2));
    for t = 1:size(Model.bodies(1).q,2)
        pos(:,t) = Model.bodies(1).q(1:3,t)+A(Model.bodies(1).q(4:7,t))*Model.bodies(1).points(i,:)';
        vel(:,t) = Model.bodies(1).qd(1:3,t)+B(Model.bodies(1).q(4:7,t),Model.bodies(1).points(i,:)')*Model.bodies(1).qd(4:7,t);
        acc(:,t) = Model.bodies(1).qdd(1:3,t)+B(Model.bodies(1).qd(4:7,t),Model.bodies(1).points(i,:)')*Model.bodies(1).qd(4:7,t)+B(Model.bodies(1).q(4:7,t),Model.bodies(1).points(i,:)')*Model.bodies(1).qdd(4:7,t);
    end

    if(i==1)
        titletext = 'Center of the Pendulum';
    else
        titletext = 'Tip of the Pendulum';
    end
    
    figure();
    subplot(2,2,1);
    plot(Model.time,pos(1,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global X Pos');
    set(gca(),'fontsize',16);
    title(titletext);   
    
    subplot(2,2,2);
    plot(Model.time,pos(2,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Y Pos');
    set(gca(),'fontsize',16);
    title(titletext);     
    
    subplot(2,2,3);
    plot(Model.time,pos(3,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Z Pos');
    set(gca(),'fontsize',16);
    title(titletext); 
    
    subplot(2,2,4);
    plot(pos(2,:),pos(3,:),'linewidth',3);
    grid on
    xlabel('Global Y Pos');
    ylabel('Global Z Pos');
    set(gca(),'fontsize',16);
    axis image
    title(titletext); 
    
    
    
    figure();
    subplot(3,1,1);
    plot(Model.time,vel(1,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global X Vel');
    set(gca(),'fontsize',16);
    title(titletext);       
    
    subplot(3,1,2);
    plot(Model.time,vel(2,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Y Vel');
    set(gca(),'fontsize',16);
    title(titletext);     
    
    subplot(3,1,3);
    plot(Model.time,vel(3,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Z Vel');
    set(gca(),'fontsize',16);
    title(titletext); 
    
    
    figure();
    subplot(3,1,1);
    plot(Model.time,acc(1,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global X Acc');
    set(gca(),'fontsize',16);
    title(titletext);       
    
    subplot(3,1,2);
    plot(Model.time,acc(2,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Y Acc');
    set(gca(),'fontsize',16);
    title(titletext);     
    
    subplot(3,1,3);
    plot(Model.time,acc(3,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Z Acc');
    set(gca(),'fontsize',16);
    title(titletext);     
end
end

function A_Matrix = A(p)
    e0 = p(1);
    e = p(2:4);
    A_Matrix = (e0^2-e'*e)*eye(3)+2*(e*e')+2*e0*Tilde(e);
end
function Cross_Product_Matrix = Tilde(v)
    Cross_Product_Matrix = [    0, -v(3),  v(2);...
                             v(3),     0, -v(1);...
                            -v(2),  v(1),     0];
end 
function B_Matrix = B(p,s)
% syms e0 e1 e2 e3 s1 s2 s3;
% p = [e0;e1;e2;e3];
% s = [s1;s2;s3];
% assume(p,'real');
% assume(s,'real');
% 
% B_Matrix = [diff(A(p)*s,e0),diff(A(p)*s,e1),diff(A(p)*s,e2),diff(A(p)*s,e3)]

e0 = p(1);
e1 = p(2);
e2 = p(3);
e3 = p(4);
s1 = s(1);
s2 = s(2);
s3 = s(3);

B_Matrix = [ 2*e0*s1 + 2*e2*s3 - 2*e3*s2, 2*e1*s1 + 2*e2*s2 + 2*e3*s3, 2*e0*s3 + 2*e1*s2 - 2*e2*s1, 2*e1*s3 - 2*e0*s2 - 2*e3*s1;...
             2*e0*s2 - 2*e1*s3 + 2*e3*s1, 2*e2*s1 - 2*e1*s2 - 2*e0*s3, 2*e1*s1 + 2*e2*s2 + 2*e3*s3, 2*e0*s1 + 2*e2*s3 - 2*e3*s2;...
             2*e0*s3 + 2*e1*s2 - 2*e2*s1, 2*e0*s2 - 2*e1*s3 + 2*e3*s1, 2*e3*s2 - 2*e2*s3 - 2*e0*s1, 2*e1*s1 + 2*e2*s2 + 2*e3*s3];

end