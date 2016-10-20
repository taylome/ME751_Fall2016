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
    coords = zeros(3,size(Model.bodies(1).q,2));
    
    for t = 1:size(Model.bodies(1).q,2)
        coords(:,t) = Model.bodies(1).q(1:3,t)+A(Model.bodies(1).q(4:7,t))*Model.bodies(1).points(i,:)';
    end

    if(i==1)
        titletext = 'Center of the Pendulum';
    else
        titletext = 'Tip of the Pendulum';
    end
    
    figure();
    subplot(2,2,1);
    plot(Model.time,coords(1,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global X Pos');
    set(gca(),'fontsize',16);
    title(titletext);   
    
    subplot(2,2,2);
    plot(Model.time,coords(2,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Y Pos');
    set(gca(),'fontsize',16);
    title(titletext);     
    
    subplot(2,2,3);
    plot(Model.time,coords(3,:),'linewidth',3);
    grid on
    xlabel('Time (s)');
    ylabel('Global Z Pos');
    set(gca(),'fontsize',16);
    title(titletext); 
    
    subplot(2,2,4);
    plot(coords(2,:),coords(3,:),'linewidth',3);
    grid on
    xlabel('Global Y Pos');
    ylabel('Global Z Pos');
    set(gca(),'fontsize',16);
    axis image
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