classdef Constraint_P
    %Implementation of the Euler Parameter Constraint
    %ME751 - Homework #6 - Oct 2016
    properties(Constant)
        type = 'P';
        numbodies = 1;        
    end
    properties
        name = '';
        id = [];
        body1 = [];
    end
    
    methods(Access = public)
        %Constructor
        function obj = Constraint_P(json_constraint)
                       
            obj.name = json_constraint.name;
            obj.id = json_constraint.id;
            obj.body1 = json_constraint.body1;

        end
        
        function Phi_P = Phi(obj, ~, q, ~)
            
            %ri = q(1:3);
            pi = q(4:7);
            
            Phi_P = (pi'*pi)-1;
            
        end
        
        function Phi_qri_P = Phi_qri(obj, ~, ~, ~)
            Phi_qri_P = zeros(1,3);
        end     
        function Phi_qpi_P = Phi_qpi(obj, ~, q, ~)
            %Phi_qpi_P = 2*pi_T

            %ri = q(1:3);
            pi = q(4:7);

            Phi_qpi_P = 2*pi';
        end
        function Phi_qi_P = Phi_qi(obj, t, q, qd)
            Phi_qi_P = [obj.Phi_qri(t, q, qd),obj.Phi_qpi(t, q, qd)];
        end
        
        function Out = Phi_qri_lambda_qri(obj,q,lambda)
            Out = zeros(3);
        end
        function Out = Phi_qri_lambda_qpi(obj,q,lambda)
            Out = zeros(3,4);
        end
        function Out = Phi_qpi_lambda_qri(obj,q,lambda)
            Out = zeros(4,3);
        end
        function Out = Phi_qpi_lambda_qpi(obj,q,lambda)
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Out = lambda*2*eye(4);
        end

        function Out = Phi_qi_lambda_qi(obj,q,lambda)
            Out = [Phi_qri_lambda_qri(obj,q,lambda),Phi_qri_lambda_qpi(obj,q,lambda);...
                   Phi_qpi_lambda_qri(obj,q,lambda),Phi_qpi_lambda_qpi(obj,q,lambda)];
        end       
        
        function Nu_P = Nu(obj, ~, ~, ~)
            Nu_P = 0;    
        end
        function Gamma_P = Gamma(obj, ~, ~, qd)
            
            %ri_d = qd(1:3);
            pi_d = qd(4:7);
            
            Gamma_P = -2*(pi_d'*pi_d);
            
        end        
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
function G_Matrix = G(p)
    e0 = p(1);
    e = p(2:4);
    G_Matrix = [-e,(-Tilde(e)+e0*eye(3))];
end
function T_Matrix = T(a)
    T_Matrix = [0,-a';...
                a,-Tilde(a)];
end
function K_Matrix = K(a,b)
    K_Matrix = 2*[a'*b, a'*Tilde(b);...
                  Tilde(a)*b, a*b'+b*a'-a'*b*eye(3)];
end
