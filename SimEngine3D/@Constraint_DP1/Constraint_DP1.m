classdef Constraint_DP1
    %Implementation of the DP1 Constraint
    %ME751 - Homework #5 - Oct 2016
    properties(Constant)
        type = 'DP1';
        numbodies = 2;        
    end
    properties
        name = '';
        id = [];
        body1 = [];
        aP1 = [];
        body2 = [];
        aP2 = [];
        func = @(t)0;
    end
    properties(Access = private)
        funcd = @(t)0;
        funcdd = @(t)0;
        has_fun = true;
    end
    
    methods(Access = public)
        %Constructor
        function obj = Constraint_DP1(json_constraint)
            
            %ensure aP1 & aP2 are columns, not rows & have a magnitude of 1
            if(size(json_constraint.aP1,1)==1)
                json_constraint.aP1 = json_constraint.aP1';
            end
            if(size(json_constraint.aP2,1)==1)
                json_constraint.aP2 = json_constraint.aP2';
            end
            
            json_constraint.aP1 = json_constraint.aP1/norm(json_constraint.aP1);
            json_constraint.aP2 = json_constraint.aP2/norm(json_constraint.aP2);
            
            
            obj.name = json_constraint.name;
            obj.id = json_constraint.id;
            obj.body1 = json_constraint.body1;
            obj.aP1 = json_constraint.aP1;
            obj.body2 = json_constraint.body2;
            obj.aP2 = json_constraint.aP2;

            
            if(strcmpi(json_constraint.fun,'NONE'))
                obj.has_fun = false;
            else
                syms t;
                symfunc = sym(eval(json_constraint.fun));
                obj.func = matlabFunction(symfunc,'vars',t);
                obj.funcd = matlabFunction(diff(symfunc),'vars',t);
                obj.funcdd = matlabFunction(diff(symfunc,2),'vars',t);
            end
        end
        
        function Phi_DP1 = Phi(obj, t, q, ~)
            %Phi_DP1 = ai_bar_T*Ai_T*Aj*aj_bar_T -f(t)
            
            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_DP1 = obj.aP1'*A(pi)'*A(pj)*obj.aP2;
            
            if(obj.has_fun) %function to define a driver
                Phi_DP1 = Phi_DP1-obj.func(t);
            end
        end
        
        function Phi_qri_DP1 = Phi_qri(obj, ~, ~, ~)
            %Phi_qri_DP1 = 0
            
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Phi_qri_DP1 = zeros(1,3);
        end
        function Phi_qrj_DP1 = Phi_qrj(obj, ~, ~, ~)
            %Phi_qrj_DP1 = 0
            
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Phi_qrj_DP1 = zeros(1,3);
        end        
        function Phi_qpi_DP1 = Phi_qpi(obj, ~, q, ~)
            %Phi_qpi_DP1 = aj_bar_T*Aj_T*B(pi,ai_bar)

            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qpi_DP1 = obj.aP2'*A(pj)'*B(pi,obj.aP1);
        end        
        function Phi_qpj_DP1 = Phi_qpj(obj, ~, q, ~)
            %Phi_qpj_DP1 = ai_bar_T*Ai_T*B(pj,aj_bar)

            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qpj_DP1 = obj.aP1'*A(pi)'*B(pj,obj.aP2);
        end
        function Phi_qi_PD1 = Phi_qi(obj, t, q, qd)
            Phi_qi_PD1 = [obj.Phi_qri(t, q, qd),obj.Phi_qpi(t, q, qd)];
        end
        function Phi_qj_DP1 = Phi_qj(obj, t, q, qd)
            Phi_qj_DP1 = [obj.Phi_qrj(t, q, qd),obj.Phi_qpj(t, q, qd)];
        end
        
        function Out = Phi_qri_lambda_qri(obj,q,lambda)
            Out = zeros(3);
        end
        function Out = Phi_qri_lambda_qrj(obj,q,lambda)
            Out = zeros(3);
        end
        function Out = Phi_qri_lambda_qpi(obj,q,lambda)
            Out = zeros(3,4);
        end
        function Out = Phi_qri_lambda_qpj(obj,q,lambda)
            Out = zeros(3,4);
        end
        
        function Out = Phi_qrj_lambda_qri(obj,q,lambda)
            Out = zeros(3);
        end
        function Out = Phi_qrj_lambda_qrj(obj,q,lambda)
            Out = zeros(3);
        end
        function Out = Phi_qrj_lambda_qpi(obj,q,lambda)
            Out = zeros(3,4);
        end
        function Out = Phi_qrj_lambda_qpj(obj,q,lambda)
            Out = zeros(3,4);
        end
        
        function Out = Phi_qpi_lambda_qri(obj,q,lambda)
            Out = zeros(4,3);
        end
        function Out = Phi_qpi_lambda_qrj(obj,q,lambda)
            Out = zeros(4,3);
        end
        function Out = Phi_qpi_lambda_qpi(obj,q,lambda)
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Out = lambda*K(obj.aP1,A(pj)*obj.aP2);
        end
        function Out = Phi_qpi_lambda_qpj(obj,q,lambda)
            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Out = lambda*B(pi,obj.aP1)'*B(pj,obj.aP2);
        end
        
        function Out = Phi_qpj_lambda_qri(obj,q,lambda)
            Out = zeros(4,3);
        end
        function Out = Phi_qpj_lambda_qrj(obj,q,lambda)
            Out = zeros(4,3);
        end
        function Out = Phi_qpj_lambda_qpi(obj,q,lambda)
            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Out = lambda*B(pj,obj.aP2)'*B(pi,obj.aP1);
        end
        function Out = Phi_qpj_lambda_qpj(obj,q,lambda)
            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Out = lambda*K(obj.aP2,A(pi)*obj.aP1);
        end
                
        function Nu_DP1 = Nu(obj, t, ~, ~)
            %Nu_DP1 = f_dot(t)
            
            Nu_DP1 = obj.funcd(t);    
        end
        function Gamma_DP1 = Gamma(obj, t, q, qd)
            %Gamma_DP1 = -2*pj_dot_T*B(pj,a_j_bar)_T*B(pi,a_i_bar)*pi_dot
            %   -aj_bar_T*Aj_T*B(pi_dot,ai_bar)*pi_dot
            %   -ai_bar_T*Ai_T*B(pj_dot,aj_bar)*pj_dot + f_dd(t)
            
            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            %ri_d = qd(1:3);
            pi_d = qd(4:7);
            %rj_d = qd(7+(1:3));
            pj_d = qd(7+(4:7));
            
            Gamma_DP1 = -2*pj_d'*B(pj,obj.aP2)'*B(pi,obj.aP1)*pi_d ...
                        -obj.aP2'*A(pj)'*B(pi_d,obj.aP1)*pi_d ...
                        -obj.aP1'*A(pi)'*B(pj_d,obj.aP2)*pj_d ...
                        +obj.funcdd(t);
            
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
