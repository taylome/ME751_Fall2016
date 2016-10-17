clc
clear all
close all

json_constraint.name = 'test';
json_constraint.id = 1;
json_constraint.body1 = 0;
json_constraint.marker1.P = [0;0;0];
json_constraint.marker1.Q = [1;0;0];
json_constraint.marker1.R = [0;0;-1];
json_constraint.body2 = 1;
json_constraint.marker2.P = [-1;0;0];
json_constraint.marker2.Q = [-1;0;1];
json_constraint.marker2.R = [0;0;0];
json_constraint.fun = 'pi()/4*cos(2*t)';

t = 0;
q = [0;0;0;...
     1;0;0;0;...
     0;0;-1;...
     cos((pi()/2)/2);[0;1;0]*sin((pi()/2)/2)];
qd = zeros(14,1);


Rev = Constraint_Revolute(json_constraint);
disp('Constraint_Revolute Test:');
Rev.Phi(t,q,qd)
Rev.Phi_qri(t,q,qd)
Rev.Phi_qrj(t,q,qd)
Rev.Phi_qpi(t,q,qd)
Rev.Phi_qpj(t,q,qd)
Rev.Nu(t,q,qd)
Rev.Gamma(t,q,qd)

mrk1 = Rev.marker1
mrk2 = Rev.marker2

qa.e0 = cos((pi()/2)/2);
e = [0;1;0]*sin((pi()/2)/2);
qa.e1 = e(1);
qa.e2 = e(2);
qa.e3 = e(3);

qb.e0 = cos((pi()/4)/2);
e = [0;0;1]*sin((pi()/4)/2);
qb.e1 = e(1);
qb.e2 = e(2);
qb.e3 = e(3);

q(1) = qa.e0 * qb.e0 - qa.e1 * qb.e1 - qa.e2 * qb.e2 - qa.e3 * qb.e3;
q(2) = qa.e0 * qb.e1 + qa.e1 * qb.e0 - qa.e3 * qb.e2 + qa.e2 * qb.e3;
q(3) = qa.e0 * qb.e2 + qa.e2 * qb.e0 + qa.e3 * qb.e1 - qa.e1 * qb.e3;
q(4) = qa.e0 * qb.e3 + qa.e3 * qb.e0 - qa.e2 * qb.e1 + qa.e1 * qb.e2;

r = [0;1*sin(pi()/4);-1*cos(pi()/4)];
q = [0;0;0;...
     1;0;0;0;...
     r;...
     q];
 
disp('Constraint_Revolute Test:');
Rev.Phi(t,q,qd)
Rev.Phi_qri(t,q,qd)
Rev.Phi_qrj(t,q,qd)
Rev.Phi_qpi(t,q,qd)
Rev.Phi_qpj(t,q,qd)
Rev.Nu(t,q,qd)
Rev.Gamma(t,q,qd)

mrk1 = Rev.marker1
mrk2 = Rev.marker2


