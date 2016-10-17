classdef Constraint_D
    %Implementation of the D Constraint
    %ME751 - Homework #5 - Oct 2016
    properties(Constant)
        type = 'D';
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
        function obj = Constraint_D(json_constraint)
            
            %ensure aP1 & aP2 are columns, not rows & have a magnitude of 1
            if(size(json_constraint.sP1,1)==1)
                json_constraint.sP1 = json_constraint.sP1';
            end
            if(size(json_constraint.sP2,1)==1)
                json_constraint.sP2 = json_constraint.sP2';
            end
            
            
            obj.name = json_constraint.name;
            obj.id = json_constraint.id;
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
        function Phi_D = Phi(obj, t, q, ~)
            %Phi_D = dij_T*dij -f(t)
            
            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_D = dij(ri,pi,obj.sP1,rj,pj,obj.sP2)'*dij(ri,pi,obj.sP1,rj,pj,obj.sP2);
            
            if(obj.has_fun) %function to define a driver
                Phi_D = Phi_D-obj.func(t);
            end
        end
        function Phi_qri_D = Phi_qri(obj, ~, q, ~)
            %Phi_qri_D = -2*dij_T
            
            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qri_D = -2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)';
        end
        function Phi_qrj_D = Phi_qrj(obj, ~, q, ~)
            %Phi_qrj_D = 2*dij_T
            
            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qrj_D = 2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)';
        end        
        function Phi_qpi_D = Phi_qpi(obj, ~, q, ~)
            %Phi_qpi_D = -2*dij_T*B(pi,si_bar_P)

            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qpi_D = -2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)'*B(pi,obj.sP1);
        end        
        function Phi_qpj_D = Phi_qpj(obj, ~, q, ~)
            %Phi_qpj_D = 2*dij_T*B(pj,sj_bar_Q)

            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            Phi_qpj_D = 2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)'*B(pj,obj.sP2);
        end         
        function Nu_D = Nu(obj, t, ~, ~)
            %Nu_D = f_dot(t)
            
            Nu_D = obj.funcd(t);    
        end
        function Gamma_D = Gamma(obj, t, q, qd)
            %Gamma_CD = -2dij_d_T*dij-2dij_T*B(pj_dot,sj_bar_Q)*pj_dot
            %   + 2dij_T*B(pi_dot,si_bar_P)*pi_dot
            %   + f_dd(t)
            
            ri = q(1:3);
            pi = q(4:7);
            rj = q(7+(1:3));
            pj = q(7+(4:7));
            
            ri_d = qd(1:3);
            pi_d = qd(4:7);
            rj_d = qd(7+(1:3));
            pj_d = qd(7+(4:7));
            
            Gamma_D = -2*(rj_d+B(pj,obj.sP2)*pj_d-ri_d-B(pi,obj.sP1)*pi_d)'*dij(ri,pi,obj.sP1,rj,pj,obj.sP2) ...
                      -2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)'*B(pj_d,obj.sP2)*pj_d ...
                      +2*dij(ri,pi,obj.sP1,rj,pj,obj.sP2)'*B(pi_d,obj.sP1)*pi_d ...
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
function d = dij(ri,pi,s_bar_p,rj,pj,s_bar_q)
    d = (rj + A(pj)*s_bar_q - ri -A(pi)*s_bar_p);
end