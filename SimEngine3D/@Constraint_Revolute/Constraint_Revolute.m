classdef Constraint_Revolute
    %Implementation of a Revolute Joint
    %ME751 - Homework #6 - Oct 2016
    properties(Constant)
        type = 'Revolute';
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
        has_fun = true;
        
        joint_angle = 0;
        
        Constraint_CD_1 = Constraint_CD.empty();
        Constraint_CD_2 = Constraint_CD.empty();
        Constraint_CD_3 = Constraint_CD.empty();
        Constraint_DP1_1 = Constraint_DP1.empty();
        Constraint_DP1_2 = Constraint_DP1.empty();
        Constraint_DP1_driver = Constraint_DP1.empty();
        Constraint_DP1_driver_alt = Constraint_DP1.empty();
        
    end
    
    methods(Access = public)
        %Constructor
        function obj = Constraint_Revolute(json_constraint)
                       
            
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
            
            %Setup the 2 DP1 constraints (x & y axis of marker 2 must be
            %perpendicular to the z axis of marker 1)
            sub_constraint.aP1 = obj.marker1.z_axis;
            sub_constraint.aP2 = obj.marker2.x_axis;
            obj.Constraint_DP1_1 = Constraint_DP1(sub_constraint);
            
            sub_constraint.aP2 = obj.marker2.y_axis;
            obj.Constraint_DP1_2 = Constraint_DP1(sub_constraint);
                       
            if(strcmpi(json_constraint.fun,'NONE'))
                obj.has_fun = false;
            else %A rotational driver was added to the joint                            
                %Setup 2 additonal DP1 constraints for the revolute joint free
                %rotational DOF driver
                %The primary constraint is PD1 marker 1 x and marker 2 x.
                %The secondary constraint when the dot product of two x axes is
                %close to 1 is PD1 marker 1 x and marker 2 y.
                sub_constraint.aP1 = obj.marker1.x_axis;
                sub_constraint.aP2 = obj.marker2.x_axis;
                sub_constraint.fun = ['cos(',json_constraint.fun,')'];
                obj.Constraint_DP1_driver = Constraint_DP1(sub_constraint);
                
                sub_constraint.aP2 = obj.marker2.y_axis;
                sub_constraint.fun = ['cos(',json_constraint.fun,'+pi()/2)'];
                obj.Constraint_DP1_driver_alt = Constraint_DP1(sub_constraint);
            end
        end
        function Phi_Rev = Phi(obj, t, q, qd)
            Phi_Rev = [obj.Constraint_CD_1.Phi( t, q, qd);...
                       obj.Constraint_CD_2.Phi( t, q, qd);...
                       obj.Constraint_CD_3.Phi( t, q, qd);...
                       obj.Constraint_DP1_1.Phi(t, q, qd);...
                       obj.Constraint_DP1_2.Phi(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                
                Phi_Rev = [Phi_Rev;...
                           obj.Constraint_DP1_driver.Phi(t, q, qd)];
            end
        end
        function Phi_qri_Rev = Phi_qri(obj, t, q, qd)
            Phi_qri_Rev = [obj.Constraint_CD_1.Phi_qri( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qri( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qri( t, q, qd);...
                           obj.Constraint_DP1_1.Phi_qri(t, q, qd);...
                           obj.Constraint_DP1_2.Phi_qri(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Phi_qri_Rev = [Phi_qri_Rev;...
                               obj.Constraint_DP1_driver.Phi_qri(t, q, qd)];
            end
        end
        function Phi_qrj_Rev = Phi_qrj(obj, t, q, qd)
            Phi_qrj_Rev = [obj.Constraint_CD_1.Phi_qrj( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qrj( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qrj( t, q, qd);...
                           obj.Constraint_DP1_1.Phi_qrj(t, q, qd);...
                           obj.Constraint_DP1_2.Phi_qrj(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Phi_qrj_Rev = [Phi_qrj_Rev;...
                               obj.Constraint_DP1_driver.Phi_qrj(t, q, qd)];
            end
        end        
        function Phi_qpi_Rev = Phi_qpi(obj, t, q, qd)
            Phi_qpi_Rev = [obj.Constraint_CD_1.Phi_qpi( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qpi( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qpi( t, q, qd);...
                           obj.Constraint_DP1_1.Phi_qpi(t, q, qd);...
                           obj.Constraint_DP1_2.Phi_qpi(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Phi_qpi_Rev = [Phi_qpi_Rev;...
                               obj.Constraint_DP1_driver.Phi_qpi(t, q, qd)];
                           
            end
        end        
        function Phi_qpj_Rev = Phi_qpj(obj, t, q, qd)
            Phi_qpj_Rev = [obj.Constraint_CD_1.Phi_qpj( t, q, qd);...
                           obj.Constraint_CD_2.Phi_qpj( t, q, qd);...
                           obj.Constraint_CD_3.Phi_qpj( t, q, qd);...
                           obj.Constraint_DP1_1.Phi_qpj(t, q, qd);...
                           obj.Constraint_DP1_2.Phi_qpj(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Phi_qpj_Rev = [Phi_qpj_Rev;...
                               obj.Constraint_DP1_driver.Phi_qpj(t, q, qd)];
                           
            end
        end         
        function Nu_Rev = Nu(obj, t, q, qd)
            Nu_Rev = [obj.Constraint_CD_1.Nu( t, q, qd);...
                           obj.Constraint_CD_2.Nu( t, q, qd);...
                           obj.Constraint_CD_3.Nu( t, q, qd);...
                           obj.Constraint_DP1_1.Nu(t, q, qd);...
                           obj.Constraint_DP1_2.Nu(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Nu_Rev = [Nu_Rev;...
                          obj.Constraint_DP1_driver.Nu(t, q, qd)];
                           
            end   
        end
        function Gamma_Rev = Gamma(obj, t, q, qd)
            Gamma_Rev = [obj.Constraint_CD_1.Gamma( t, q, qd);...
                           obj.Constraint_CD_2.Gamma( t, q, qd);...
                           obj.Constraint_CD_3.Gamma( t, q, qd);...
                           obj.Constraint_DP1_1.Gamma(t, q, qd);...
                           obj.Constraint_DP1_2.Gamma(t, q, qd)];
            
            if(obj.has_fun) %function to define a driver along the otherwise free rotational DOF
                Gamma_Rev = [Gamma_Rev;...
                          obj.Constraint_DP1_driver.Gamma(t, q, qd)];
                           
            end
        end
    end
    
end
