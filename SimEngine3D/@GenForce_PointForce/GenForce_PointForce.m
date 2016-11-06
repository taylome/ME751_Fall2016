classdef GenForce_PointForce
    %Implementation of the Point Force Generalized Force
    %ME751 - Homework #7 - Oct 2016
    properties(Constant)
        type = 'PointForce';
        numbodies = 1;
    end
    properties
        name = '';
        id = [];
        body1 = [];
        sP1 = [];
        GlobalFrame = true;
        funcX = @(t)0;
        funcY = @(t)0;
        funcZ = @(t)0;
    end
    
    methods
        %Constructor
        function obj = GenForce_PointForce(json_force)
            %Ensure fun ~= 'NONE'
            if(strcmpi(json_force.funX,'NONE'))
                error('Point Force Generalized Force must have a value for "funX" other than "NONE"');
            end
            if(strcmpi(json_force.funY,'NONE'))
                error('Point Force Generalized Force must have a value for "funY" other than "NONE"');
            end
            if(strcmpi(json_force.funZ,'NONE'))
                error('Point Force Generalized Force must have a value for "funY" other than "NONE"');
            end
            
            %ensure sP1 is a column, not a row
            if(size(json_force.sP1,1)==1)
                json_force.sP1 = json_force.sP1';
            end
            
            obj.name = json_force.name;
            obj.id = json_force.id;
            obj.body1 = json_force.body1;
            obj.sP1 = json_force.sP1;
            if(strcmpi(json_force.frame,'GRF'))
                obj.GlobalFrame = true;
            elseif(strcmpi(json_force.frame,'LRF'))
                obj.GlobalFrame = false;
            else
                error('Point Force Generalized Force frame must be either "GRF" or "LRF"');
            end
            
            syms t;
            obj.funcX = matlabFunction(sym(eval(json_force.funX)),'vars',t);
            obj.funcY = matlabFunction(sym(eval(json_force.funY)),'vars',t);
            obj.funcZ = matlabFunction(sym(eval(json_force.funZ)),'vars',t);
        end
        function Q_PntFrc = Q(obj, t, q, ~)
            pi = q(4:7);
            
            F = [obj.funcX(t);obj.funcY(t);obj.funcZ(t)];
            if(~obj.GlobalFrame)
                F = A(pi)*F;
            end
            
            Q_PntFrc = [F;2*G(pi)'*Tilde(obj.sP1)*A(pi)'*F];
        end
        
        function Out = Q_ri(obj,t,q,qd)
            Out = zeros(7,3);
        end
        function Out = Q_rdi(obj,t,q,qd)
            Out = zeros(7,3);
        end
        function Out = Q_pi(obj,t,q,qd)
            %Calculate Numerically for now
            delta = 0.0001;
            Out = zeros(7,4);
            QBase = Q(obj, t, q, qd);
            for i = 1:4
                qdelta = zeros(7,1);
                qdelta(4+i) = delta;
                Out(:,i) = (Q(obj, t, q+qdelta, qd)-QBase)/delta;
            end
        end        
        function Out = Q_pdi(obj,t,q,qd)
            Out = zeros(7,4);
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
