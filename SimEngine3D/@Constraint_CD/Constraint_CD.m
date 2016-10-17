classdef Constraint_CD
    %Implementation of the CD Constraint
    %ME751 - Homework #5 - Oct 2016
    properties(Constant)
        type = 'CD';
        numbodies = 2;        
    end
    properties
        name = '';
        id = [];
        c = [];
        body1 = [];
        sP1 = [];
        body2 = [];
        sP2 = [];
        func = @(t)0;
    end
    properties(Access = private)
        funcd = @(t)0;
        funcdd = @(t)0;
        has_fun = true;
    end
    
    methods(Access = public)
        %Constructor
        function obj = Constraint_CD(json_constraint)
            
            %ensure aP1 & aP2 are columns, not rows & have a magnitude of 1
            if(size(json_constraint.sP1,1)==1)
                json_constraint.sP1 = json_constraint.sP1';
            end
            if(size(json_constraint.sP2,1)==1)
                json_constraint.sP2 = json_constraint.sP2';
            end
            if(size(json_constraint.c,1)==1)
                json_constraint.c = json_constraint.c';
            end            
            json_constraint.c = json_constraint.c/norm(json_constraint.c);
            
            
            obj.name = json_constraint.name;
            obj.id = json_constraint.id;
            obj.c = json_constraint.c;
            obj.body1 = json_constraint.body1;
            obj.sP1 = json_constraint.sP1;
            obj.body2 = json_constraint.body2;
            obj.sP2 = json_constraint.sP2;

            
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
        function Phi_CD = Phi(obj, t, q, ~)
            %Phi_CD = c_T*dij -f(t)
            
            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_CD = obj.c'*(rj+A(pj)*obj.sP2-ri-A(pi)*obj.sP1);
            
            if(obj.has_fun) %function to define a driver
                Phi_CD = Phi_CD-obj.func(t);
            end
        end
        function Phi_qri_CD = Phi_qri(obj, ~, ~, ~)
            %Phi_qri_CD = -c_T
            
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Phi_qri_CD = -obj.c';
        end
        function Phi_qrj_CD = Phi_qrj(obj, ~, ~, ~)
            %Phi_qrj_CD = c_T
            
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Phi_qrj_CD = obj.c';
        end        
        function Phi_qpi_CD = Phi_qpi(obj, ~, q, ~)
            %Phi_qpi_CD = -c_T*B(pi,si_bar_P)

            %ri = q(1:3);
            pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            Phi_qpi_CD = -obj.c'*B(pi,obj.sP1);
        end        
        function Phi_qpj_CD = Phi_qpj(obj, ~, q, ~)
            %Phi_qpj_CD = c_T*B(pj,sj_bar_Q)

            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qpj_CD = obj.c'*B(pj,obj.sP2);
        end         
        function Nu_CD = Nu(obj, t, ~, ~)
            %Nu_CD = f_dot(t)
            
            Nu_CD = obj.funcd(t);    
        end
        function Gamma_CD = Gamma(obj, t, ~, qd)
            %Gamma_CD = c_T*B(pi_dot,si_bar_P)*pi_dot
            %   - c_T*B(pj_dot,sj_bar_Q)*pj_dot
            %   + f_dd(t)
            
            %ri = q(1:3);
            %pi = q(4:7);
            %rj = q(7+(1:3));
            %pj = q(7+(4:7));
            
            %ri_d = qd(1:3);
            pi_d = qd(4:7);
            %rj_d = qd(7+(1:3));
            pj_d = qd(7+(4:7));
            
            Gamma_CD = obj.c'*B(pi_d,obj.sP1)*pi_d ...
                      -obj.c'*B(pj_d,obj.sP2)*pj_d ...
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
