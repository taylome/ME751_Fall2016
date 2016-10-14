clc
clear all
close all

json_constraint.name = 'test';
json_constraint.id = 1;
json_constraint.body1 = 1;
json_constraint.aP1 = [1;0;0];
json_constraint.sP1 = [0;0;1];
json_constraint.body2 = 2;
json_constraint.aP2 = [1;0;0];
json_constraint.sP2 = [0;0;1];
json_constraint.fun = 'NONE';

t = 0;
q = [0;0;0;...
     1;0;0;0;...
     0;0;0;...
     1;0;0;0];
qd = zeros(14,1);


DP1 = Constraint_DP1(json_constraint);
disp('DP1 Test:');
DP1.Phi(t,q,qd)
DP1.Phi_qri(t,q,qd)
DP1.Phi_qrj(t,q,qd)
DP1.Phi_qpi(t,q,qd)
DP1.Phi_qpj(t,q,qd)
DP1.Nu(t,q,qd)
DP1.Gamma(t,q,qd)

DP2 = Constraint_DP2(json_constraint);
disp('DP2 Test:');
DP2.Phi(t,q,qd)
DP2.Phi_qri(t,q,qd)
DP2.Phi_qrj(t,q,qd)
DP2.Phi_qpi(t,q,qd)
DP2.Phi_qpj(t,q,qd)
DP2.Nu(t,q,qd)
DP2.Gamma(t,q,qd)

D = Constraint_D(json_constraint);
disp('D Test:');
D.Phi(t,q,qd)
D.Phi_qri(t,q,qd)
D.Phi_qrj(t,q,qd)
D.Phi_qpi(t,q,qd)
D.Phi_qpj(t,q,qd)
D.Nu(t,q,qd)
D.Gamma(t,q,qd)

CD = Constraint_D(json_constraint);
disp('CD Test:');
CD.Phi(t,q,qd)
CD.Phi_qri(t,q,qd)
CD.Phi_qrj(t,q,qd)
CD.Phi_qpi(t,q,qd)
CD.Phi_qpj(t,q,qd)
CD.Nu(t,q,qd)
CD.Gamma(t,q,qd)