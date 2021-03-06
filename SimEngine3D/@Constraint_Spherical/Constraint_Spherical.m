classdef Constraint_Spherical
    %Implementation of a Spherical Joint
    %ME751 - Homework #8 - Nov 2016
    properties(Constant)
        type = 'Spherical';
        numbodies = 2;        
    end
    properties
        name = '';
        id = [];
        body1 = [];
        marker1 = struct();
        body2 = [];
        marker2 = struct();
        func = @(t)0;
    end
    properties(Access = private)
        funcd = @(t)0;
        funcdd = @(t)0;
        has_fun = false;
        
        joint_angle = 0;
        
        Constraint_CD_1 = Constraint_CD.empty();
        Constraint_CD_2 = Constraint_CD.empty();
        Constraint_CD_3 = Constraint_CD.empty();
        
    end
    
    methods(Access = public)
        %Constructor
        function obj = Constraint_Spherical(json_constraint)
                       
            
            obj.name = json_constraint.name;
            obj.id = json_constraint.id;
            
            %P = origin
            %Q = point along the z axis
            %R = point in the XZ plane (positive x coordinate)
            obj.body1 = json_constraint.body1;
            z_axis = json_constraint.marker1.Q-json_constraint.marker1.P;
            z_axis = z_axis/norm(z_axis);
            r_axis = json_constraint.marker1.R-json_constraint.marker1.P;
            r_axis = r_axis/norm(r_axis);
            y_axis = cross(z_axis,r_axis);
            y_axis = y_axis/norm(y_axis);
            x_axis = cross(y_axis,z_axis);
            x_axis = x_axis/norm(x_axis);
            obj.marker1.origin = json_constraint.marker1.P;
            obj.marker1.x_axis = x_axis;
            obj.marker1.y_axis = y_axis;
            obj.marker1.z_axis = z_axis;
            
            obj.body2 = json_constraint.body2;
            z_axis = json_constraint.marker2.Q-json_constraint.marker2.P;
            z_axis = z_axis/norm(z_axis);
            r_axis = json_constraint.marker2.R-json_constraint.marker2.P;
            r_axis = r_axis/norm(r_axis);
            y_axis = cross(z_axis,r_axis);
            y_axis = y_axis/norm(y_axis);
            x_axis = cross(y_axis,z_axis);
            x_axis = x_axis/norm(x_axis);
            obj.marker2.origin = json_constraint.marker2.P;
            obj.marker2.x_axis = x_axis;
            obj.marker2.y_axis = y_axis;
            obj.marker2.z_axis = z_axis;
            
            
            
            %Set up the 3 CD joints
            sub_constraint = json_constraint;
            sub_constraint.sP1 = json_constraint.marker1.P;
            sub_constraint.sP2 = json_constraint.marker2.P;
            sub_constraint.fun = 'NONE';
            sub_constraint.c = [1;0;0];
            obj.Constraint_CD_1 = Constraint_CD(sub_constraint);
            sub_constraint.c = [0;1;0];
            obj.Constraint_CD_2 = Constraint_CD(sub_constraint);
            sub_constraint.c = [0;0;1];
            obj.Constraint_CD_3 = Constraint_CD(sub_constraint);
                       
            if(strcmpi(json_constraint.fun,'NONE'))
                obj.has_fun = false;
            else %A rotational driver was added to the joint                            
                %Add in driver functionality later
                obj.has_fun = false;
            end
        end
        
        function Phi_Rev = Phi(obj, t, q, qd)
            Phi_Rev = [obj.Constraint_CD_1.Phi( t, q, qd);...
                       obj.Constraint_CD_2.Phi( t, q, qd);...
                       obj.Constraint_CD_3.Phi( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF               
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
            end
        end
        
        function Phi_qri_Rev = Phi_qri(obj, t, q, qd)
            Phi_qri_Rev = [obj.Constraint_CD_1.Phi_qri( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qri( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qri( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
            end
        end
        function Phi_qrj_Rev = Phi_qrj(obj, t, q, qd)
            Phi_qrj_Rev = [obj.Constraint_CD_1.Phi_qrj( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qrj( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qrj( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
            end
        end        
        function Phi_qpi_Rev = Phi_qpi(obj, t, q, qd)
            Phi_qpi_Rev = [obj.Constraint_CD_1.Phi_qpi( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qpi( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qpi( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
                           
            end
        end        
        function Phi_qpj_Rev = Phi_qpj(obj, t, q, qd)
            Phi_qpj_Rev = [obj.Constraint_CD_1.Phi_qpj( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qpj( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qpj( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
                           
            end
        end
        function Phi_qi_Rev = Phi_qi(obj, t, q, qd)
            Phi_qi_Rev = [obj.Phi_qri(t, q, qd),obj.Phi_qpi(t, q, qd)];
        end
        function Phi_qj_Rev = Phi_qj(obj, t, q, qd)
            Phi_qj_Rev = [obj.Phi_qrj(t, q, qd),obj.Phi_qpj(t, q, qd)];
        end  
        
        function Out = Phi_qri_lambda_qri(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qri_lambda_qri( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qri_lambda_qri( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qri_lambda_qri( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qri_lambda_qrj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qri_lambda_qrj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qri_lambda_qrj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qri_lambda_qrj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qri_lambda_qpi(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qri_lambda_qpi( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qri_lambda_qpi( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qri_lambda_qpi( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qri_lambda_qpj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qri_lambda_qpj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qri_lambda_qpj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qri_lambda_qpj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end

        function Out = Phi_qrj_lambda_qri(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qrj_lambda_qri( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qrj_lambda_qri( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qrj_lambda_qri( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qrj_lambda_qrj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qrj_lambda_qrj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qrj_lambda_qrj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qrj_lambda_qrj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qrj_lambda_qpi(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qrj_lambda_qpi( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qrj_lambda_qpi( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qrj_lambda_qpi( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qrj_lambda_qpj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qrj_lambda_qpj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qrj_lambda_qpj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qrj_lambda_qpj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        
        function Out = Phi_qpi_lambda_qri(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpi_lambda_qri( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpi_lambda_qri( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpi_lambda_qri( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qpi_lambda_qrj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpi_lambda_qrj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpi_lambda_qrj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpi_lambda_qrj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qpi_lambda_qpi(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpi_lambda_qpi( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpi_lambda_qpi( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpi_lambda_qpi( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        function Out = Phi_qpi_lambda_qpj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpi_lambda_qpj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpi_lambda_qpj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpi_lambda_qpj( q,lambda(3));
            if(obj.has_fun)
                  
            end
        end
        
        function Out = Phi_qpj_lambda_qri(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpj_lambda_qri( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpj_lambda_qri( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpj_lambda_qri( q,lambda(3));
            if(obj.has_fun)

            end
        end
        function Out = Phi_qpj_lambda_qrj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpj_lambda_qrj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpj_lambda_qrj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpj_lambda_qrj( q,lambda(3));
            if(obj.has_fun)
                
            end
        end
        function Out = Phi_qpj_lambda_qpi(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpj_lambda_qpi( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpj_lambda_qpi( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpj_lambda_qpi( q,lambda(3));
            if(obj.has_fun)
                
            end
        end
        function Out = Phi_qpj_lambda_qpj(obj,q,lambda)
            Out = obj.Constraint_CD_1.Phi_qpj_lambda_qpj( q,lambda(1)) + ...
                  obj.Constraint_CD_2.Phi_qpj_lambda_qpj( q,lambda(2)) + ...
                  obj.Constraint_CD_3.Phi_qpj_lambda_qpj( q,lambda(3));
            if(obj.has_fun)

            end
        end

        function Out = Phi_qi_lambda_qi(obj,q,lambda)
            Out = [obj.Phi_qri_lambda_qri(q,lambda),obj.Phi_qri_lambda_qpi(q,lambda);...
                   obj.Phi_qpi_lambda_qri(q,lambda),obj.Phi_qpi_lambda_qpi(q,lambda)];
        end
        function Out = Phi_qi_lambda_qj(obj,q,lambda)
            Out = [obj.Phi_qri_lambda_qrj(q,lambda),obj.Phi_qri_lambda_qpj(q,lambda);...
                   obj.Phi_qpi_lambda_qrj(q,lambda),obj.Phi_qpi_lambda_qpj(q,lambda)];
        end
        function Out = Phi_qj_lambda_qi(obj,q,lambda)
            Out = [obj.Phi_qrj_lambda_qri(q,lambda),obj.Phi_qrj_lambda_qpi(q,lambda);...
                   obj.Phi_qpj_lambda_qri(q,lambda),obj.Phi_qpj_lambda_qpi(q,lambda)];
        end
        function Out = Phi_qj_lambda_qj(obj,q,lambda)
            Out = [obj.Phi_qrj_lambda_qrj(q,lambda),obj.Phi_qrj_lambda_qpj(q,lambda);...
                   obj.Phi_qpj_lambda_qrj(q,lambda),obj.Phi_qpj_lambda_qpj(q,lambda)];
        end
        
        function Nu_Rev = Nu(obj, t, q, qd)
            Nu_Rev = [obj.Constraint_CD_1.Nu( t, q, qd);...
                           obj.Constraint_CD_2.Nu( t, q, qd);...
                           obj.Constraint_CD_3.Nu( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
                           
            end   
        end
        function Gamma_Rev = Gamma(obj, t, q, qd)
            Gamma_Rev = [obj.Constraint_CD_1.Gamma( t, q, qd);...
                           obj.Constraint_CD_2.Gamma( t, q, qd);...
                           obj.Constraint_CD_3.Gamma( t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                %ri = q(1:3);
                pi = q(4:7);
                %rj = q(7+(1:3));
                pj = q(7+(4:7));
           
                           
            end
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

