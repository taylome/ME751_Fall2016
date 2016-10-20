%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Oct 2016
%
% Homework 6 - Problem #2
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------


clc
clear all
close all

[pathstr,~,~] = fileparts(mfilename('fullpath'));
model_name = [pathstr,'\models\revJoint.mdl'];

Model = simEngine3D(model_name);

t = 0;
q = [0;0;0;1;0;0;0;Model.bodies.q(:,1)];
qd = [zeros(7,1);Model.bodies.qd(:,1)];

PHI_Full = [Model.constraints{1}.Phi(t,q,qd);Model.constraints{2}.Phi(t,q,qd)]

PHI_q_Full = [Model.constraints{1}.Phi_qj(t,q,qd);Model.constraints{2}.Phi_qi(t,q,qd)]

NU = [Model.constraints{1}.Nu(t,q,qd);Model.constraints{2}.Nu(t,q,qd)]

GAMMA = [Model.constraints{1}.Gamma(t,q,qd);Model.constraints{2}.Gamma(t,q,qd)]

%--------------------------------------------------------------------------
%OUTPUT:
%--------------------------------------------------------------------------
% PHI_Full =
%    1.0e-15 *
%     0.0278
%          0
%     0.1110
%    -0.0278
%          0
%     0.1110
%          0
% PHI_q_Full =
%     1.0000         0         0   -1.3066   -0.5412    1.3066    0.5412
%          0    1.0000         0   -0.5412   -1.3066   -0.5412   -1.3066
%          0         0    1.0000    1.3066   -0.5412    1.3066   -0.5412
%          0         0         0    1.3066    0.5412   -1.3066   -0.5412
%          0         0         0   -0.5412    1.3066    0.5412   -1.3066
%          0         0         0    1.3066   -0.5412    1.3066   -0.5412
%          0         0         0    2.0000         0         0         0
% NU =
%      0
%      0
%      0
%      0
%      0
%      0
%      0
% GAMMA =
%          0
%          0
%          0
%          0
%          0
%     2.2214
%          0