function Model = simEngine3D(model_name)
%--------------------------------------------------------------------------
% simEngine3D - ME 751
% Mike Taylor - Oct 2016
%
% Model = simEngine3D(model_name)
% model_name = string with the name of the mdl file
%
% simEngine3D REQUIRES MATLAB R2016b OR LATER
%--------------------------------------------------------------------------

tic
if(nargin<1)
    [pathstr,~,~] = fileparts(mfilename('fullpath'));
    model_name = [pathstr,'\models\me751_HW07.mdl'];
end

%--------------------------------------------------------------------------
% Settings
%--------------------------------------------------------------------------
Settings.NewtonRaphsonMaxIterations = 50;
Settings.NewtonRaphsonTol = 1E-8;
Settings.model_name = model_name;
%--------------------------------------------------------------------------
 
disp(['Reading Model: ',model_name]);

[Model] = read_Model(model_name);
%assignin('base','Model',Model)

if(strcmpi(Model.simulation.Type,'Kinematics'))
    disp('Executing Kinematic Analysis...');
    Model = Run_KinematicAnalysis(Model,Settings,false);
elseif(strcmpi(Model.simulation.Type,'Inverse Dynamics'))
    disp('Executing Inverse Dynamics Analysis...');
    Model = Run_KinematicAnalysis(Model,Settings,true);    
elseif(strcmpi(Model.simulation.Type,'Dynamics'))
    error('Dynamic Analysis has not been written yet');
end
toc()

assignin('base','Model',Model);

end

function Model = Run_KinematicAnalysis(Model,Settings,CalcReactions)

NumTimeSteps = Model.simulation.tend/Model.simulation.stepSize+1;

for i = 1:length(Model.bodies)
     Model.bodies(i).q = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
     Model.bodies(i).qd = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
     Model.bodies(i).qdd = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
end

if(CalcReactions)
    for i = 1:length(Model.constraints)
        Gamma_c = Model.constraints{i}.Gamma(0, zeros(14,1), zeros(14,1));
        Model.ConstraintReactions(i).Values = zeros(length(Gamma_c),NumTimeSteps);
    end
end

%--------------------------------------------------------------------------
%Step 1: Generate inital q vector
%Step 2: Postion Analysis - Newton-Raphson
%        A) Generate Phi & Phi_q [Jacobian]
%        B) Calculate Correction
%        C) Apply Correction
%        D) Break if Correction is less than tolerance or max iterations is
%           reached
%Step 3: Velocity Analysis - Solve Phi_q*qdot = Nu
%Step 4: Acceleration Analysis - Solve Phi_q*qdot = Gamma
%Step 5: Store the q vectors into the bodies & use for the next inital 
%        guess - Return to Step 2 until final time is reached
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Step 1: Generate inital q vector
%--------------------------------------------------------------------------
q = zeros(sum(Model.NumGeneralizedCoordinates),1);
for i = 1:length(Model.bodies)
    q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i)) = Model.bodies(i).q0;
end

for tindex = 1:NumTimeSteps
    t = (tindex-1)*Model.simulation.stepSize;
    
    if(mod(tindex,1000)==1)
        disp(['t = ',num2str(t),'s']);
    end

%--------------------------------------------------------------------------
%Step 2: Postion Analysis - Newton-Raphson
%--------------------------------------------------------------------------    
    for NR = 1:Settings.NewtonRaphsonMaxIterations 
    % A) Generate Phi & Phi_q [Jacobian]
        Phi = zeros(size(q));
        Phi_q = zeros(size(q,1));
        index = 1;

        for i = 1:length(Model.constraints)
            if(Model.constraints{i}.numbodies == 1)
                StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                qi = q(StartIndex_i:EndIndex_i);

                Phi_c = Model.constraints{i}.Phi(t, qi, []);
                
                Phi(index:index+(length(Phi_c)-1)) = Phi_c;
                Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);

                index = index + length(Phi_c);
            elseif(Model.constraints{i}.body1 == 0) %Connected to ground
                StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                qi = [[0;0;0;1;0;0;0];q(StartIndex_j:EndIndex_j)];
                
                Phi_c = Model.constraints{i}.Phi(t, qi, []);

                Phi(index:index+(length(Phi_c)-1)) = Phi_c;
                Phi_q(index:index+(length(Phi_c)-1),StartIndex_j:EndIndex_j) = Model.constraints{i}.Phi_qj(t, qi, []);
                
                index = index + length(Phi_c);
            elseif(Model.constraints{i}.body2 == 0) %Connected to ground
                StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                qi = [q(StartIndex_i:EndIndex_i);[0;0;0;1;0;0;0]];
                
                Phi_c = Model.constraints{i}.Phi(t, qi, []);
                
                Phi(index:index+(length(Phi_c)-1)) = Phi_c;
                Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);
                
                index = index + length(Phi_c);                
            else
                StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];
                
                Phi_c = Model.constraints{i}.Phi(t, qi, []);
                
                Phi(index:index+(length(Phi_c)-1)) = Phi_c;
                Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);
                Phi_q(index:index+(length(Phi_c)-1),StartIndex_j:EndIndex_j) = Model.constraints{i}.Phi_qj(t, qi, []);
                
                index = index + length(Phi_c);
            end
        end
        
        % B) Calculate Correction
        correction = Phi_q\Phi;
        
        % C) Apply Correction        
        q = q - correction;
        
        % D) Break if Correction is less than tolerance or max iterations
        %    is reached
        if(norm(correction)<Settings.NewtonRaphsonTol)
            break  
        end
    end

%--------------------------------------------------------------------------
%Step 3: Velocity Analysis - Solve Phi_q*qdot = Nu
%--------------------------------------------------------------------------
    index = 1;
    Nu = zeros(size(q));
    for i = 1:length(Model.constraints)
        if(Model.constraints{i}.numbodies == 1)
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            qi = q(StartIndex_i:EndIndex_i);
        elseif(Model.constraints{i}.body1 == 0) %Connected to ground
            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
            qi = [[0;0;0;1;0;0;0];q(StartIndex_j:EndIndex_j)];
        elseif(Model.constraints{i}.body2 == 0) %Connected to ground
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            qi = [q(StartIndex_i:EndIndex_i);[0;0;0;1;0;0;0]];        
        else
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
            qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];
        end

        Nu_c = Model.constraints{i}.Nu(t, qi, []);
        Nu(index:index+(length(Nu_c)-1)) = Nu_c;

        index = index + length(Nu_c);
    end
    qd = Phi_q\Nu;

%--------------------------------------------------------------------------
%Step 4: Acceleration Analysis - Solve Phi_q*qdot = Gamma
%--------------------------------------------------------------------------
    index = 1;
    Gamma = zeros(size(q));
    for i = 1:length(Model.constraints)
        if(Model.constraints{i}.numbodies == 1)
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            qi = q(StartIndex_i:EndIndex_i);
            qdi = qd(StartIndex_i:EndIndex_i);
        elseif(Model.constraints{i}.body1 == 0) %Connected to ground
            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
            qi = [[0;0;0;1;0;0;0];q(StartIndex_j:EndIndex_j)];
            qdi = [zeros(7,1);qd(StartIndex_j:EndIndex_j)];
        elseif(Model.constraints{i}.body2 == 0) %Connected to ground
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            qi = [q(StartIndex_i:EndIndex_i);[0;0;0;1;0;0;0]];             
            qdi = [qd(StartIndex_i:EndIndex_i);zeros(7,1)];
        else
            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
            qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];
            qdi = [qd(StartIndex_i:EndIndex_i);qd(StartIndex_j:EndIndex_j)];
        end

        Gamma_c = Model.constraints{i}.Gamma(t, qi, qdi);
        Gamma(index:index+(length(Gamma_c)-1)) = Gamma_c;

        index = index + length(Gamma_c);
    end
    qdd = Phi_q\Gamma;

%--------------------------------------------------------------------------
%Step 5: Store the q vectors into the bodies & use for the next inital 
%        guess
%--------------------------------------------------------------------------
    for i = 1:length(Model.bodies)
        Model.bodies(i).q(:,tindex) = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        Model.bodies(i).qd(:,tindex) = qd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        Model.bodies(i).qdd(:,tindex) = qdd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
    end
    

%--------------------------------------------------------------------------
%Step 6: Calculate the Inverse Dynamics Reaction Forces
%        - Return to Step 2 until final time is reached
%--------------------------------------------------------------------------
    
    if(CalcReactions)
        RHS = zeros(size(Phi_q,2),1);

        for i = 1:length(Model.bodies)
            StartIndex_i = Model.bodies(i).StartIndex;

            RHS(StartIndex_i:StartIndex_i+2) = -Model.bodies(i).mass*eye(3) * Model.bodies(i).qdd(1:3,tindex);
            RHS(StartIndex_i+3:StartIndex_i+6) = -4*G(Model.bodies(i).q(4:7,tindex))'*diag(Model.bodies(i).MOI)*G(Model.bodies(i).q(4:7,tindex)) * Model.bodies(i).qdd(4:7,tindex) + ...
                                                  8*G(Model.bodies(i).qd(4:7,tindex))'*diag(Model.bodies(i).MOI)*G(Model.bodies(i).qd(4:7,tindex)) * Model.bodies(i).q(4:7,tindex);
        end
        
        for i = 1:length(Model.forces)
            StartIndex_i = Model.bodies(Model.forces{i}.body1).StartIndex;

            RHS(StartIndex_i:StartIndex_i+6) = RHS(StartIndex_i:StartIndex_i+6) + Model.forces{i}.Q(t,Model.bodies(i).q(:,tindex),Model.bodies(i).qd(:,tindex));
        end
        
        % Solve for the Lagrange multipliers
        lambdas = Phi_q'\RHS;
        
        % Solve for the reaction forces and torques for each constraint
        index = 1;
        for i = 1:length(Model.constraints)
            AccumForceTrq_Bodyi = zeros(6,1);
            AccumForceTrq_Bodyj = zeros(6,1);
            
            NumEquations = length(Model.constraints{i}.Gamma(0, zeros(14,1), zeros(14,1)));
            
            if(Model.constraints{i}.numbodies == 1)
                %Need to add markers to all of the GCONs and P constraints
%                 StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
%                 EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
%                 qi = q(StartIndex_i:EndIndex_i);
%                 Phi_qi = Phi_q(index:(index-1+NumEquations),StartIndex_i:EndIndex_i);
%                 
%                 for c = 1:NumEquations
%                     lambda = lambdas(index + c-1);
%                     Force = -Phi_qi(c,1:3)'*lambda; %Global Frame
%                     Torque = -1/2*G(qi(4:7))*Phi_qi(c,4:7)'*lambda-Tilde(Model.constraints{i}.marker1.origin)*A(qi(4:7))'*Force;
%                     Torque = A(qi(4:7))*Torque; %Local -> Global Frame
%                 end
                
            elseif(Model.constraints{i}.body1 == 0) %Connected to ground
                StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                qj = q(StartIndex_j:EndIndex_j);
                Phi_qj = Phi_q(index:(index-1+NumEquations),StartIndex_j:EndIndex_j);
                
                for c = 1:NumEquations
                    lambda = lambdas(index + c-1);
                    Force = -Phi_qj(c,1:3)'*lambda; %Global Frame
                    Torque = -1/2*G(qj(4:7))*Phi_qj(c,4:7)'*lambda-Tilde(Model.constraints{i}.marker2.origin)*A(qj(4:7))'*Force;
                    Torque = A(qj(4:7))*Torque; %Local -> Global Frame
                    
                    AccumForceTrq_Bodyj = AccumForceTrq_Bodyj + [Force;Torque];
                end
                AccumForceTrq_Bodyi = -AccumForceTrq_Bodyj;                
                
            elseif(Model.constraints{i}.body2 == 0) %Connected to ground
                StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                qi = q(StartIndex_i:EndIndex_i);
                Phi_qi = Phi_q(index:(index-1+NumEquations),StartIndex_i:EndIndex_i);
                
                for c = 1:NumEquations
                    lambda = lambdas(index + c-1);
                    Force = -Phi_qi(c,1:3)'*lambda; %Global Frame
                    Torque = -1/2*G(qi(4:7))*Phi_qi(c,4:7)'*lambda-Tilde(Model.constraints{i}.marker1.origin)*A(qi(4:7))'*Force;
                    Torque = A(qi(4:7))*Torque; %Local -> Global Frame
                    
                    AccumForceTrq_Bodyi = AccumForceTrq_Bodyi + [Force;Torque];
                end
                AccumForceTrq_Bodyj = -AccumForceTrq_Bodyi;
                
            else
                StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                qi = q(StartIndex_i:EndIndex_i);
                Phi_qi = Phi_q(index:(index-1+NumEquations),StartIndex_i:EndIndex_i);
                qj = q(StartIndex_j:EndIndex_j);
                Phi_qj = Phi_q(index:(index-1+NumEquations),StartIndex_j:EndIndex_j);
                
                for c = 1:NumEquations
                    lambda = lambdas(index + c-1);
                    Force = -Phi_qi(c,1:3)'*lambda; %Global Frame
                    Torque = -1/2*G(qi(4:7))*Phi_qi(c,4:7)'*lambda-Tilde(Model.constraints{i}.marker1.origin)*A(qi(4:7))'*Force;
                    Torque = A(qi(4:7))*Torque; %Local -> Global Frame
                    
                    AccumForceTrq_Bodyi = AccumForceTrq_Bodyi + [Force;Torque];
                    
                    Force = -Phi_qj(c,1:3)'*lambda; %Global Frame
                    Torque = -1/2*G(qj(4:7))*Phi_qj(c,4:7)'*lambda-Tilde(Model.constraints{i}.marker2.origin)*A(qj(4:7))'*Force;
                    Torque = A(qj(4:7))*Torque; %Local -> Global Frame
                    
                    AccumForceTrq_Bodyj = AccumForceTrq_Bodyj + [Force;Torque];                    
                end                

            end            

            Model.ConstraintReactions(i).Body1(:,tindex) = AccumForceTrq_Bodyi;
            Model.ConstraintReactions(i).Body2(:,tindex) = AccumForceTrq_Bodyj;
            
            index = index + NumEquations;
        end
    end
    
end

Model.time = ((0:NumTimeSteps-1)*Model.simulation.stepSize)';

% %Write visualization results file
% results = zeros(NumTimeSteps,size(q,1)+1);
% results(:,1) = Model.time;
% [~,Index] = sort(Model.BodyIDs);
% pos = 2;
% for i = 1:length(Index)
%     results(:,pos:(pos-1+size(Model.bodies(i).q,1))) = Model.bodies(i).q';
%     pos = pos+size(Model.bodies(i).q,1);
% end
% dlmwrite([Settings.model_name,'.res'],results,'delimiter',',','precision',4)
% 
% visualize([Settings.model_name,'.adm'], [Settings.model_name,'.res']);




end

function Model = Run_InverseDynamicsAnalysis(Model,Settings)

t = 0;
q = zeros(sum(Model.NumGeneralizedCoordinates),1);
for i = 1:length(Model.bodies)
    q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i)) = Model.bodies(i).q0;
end

Phi = zeros(size(q));
Phi_q = zeros(size(q,1));
index = 1;

for i = 1:length(Model.constraints)
    if(Model.constraints{i}.numbodies == 1)
        StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
        EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
        qi = q(StartIndex_i:EndIndex_i);

        Phi_c = Model.constraints{i}.Phi(t, qi, []);

        Phi(index:index+(length(Phi_c)-1)) = Phi_c;
        Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);

        index = index + length(Phi_c);
    elseif(Model.constraints{i}.body1 == 0) %Connected to ground
        StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
        EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
        qi = [[0;0;0;1;0;0;0];q(StartIndex_j:EndIndex_j)];

        Phi_c = Model.constraints{i}.Phi(t, qi, []);

        Phi(index:index+(length(Phi_c)-1)) = Phi_c;
        Phi_q(index:index+(length(Phi_c)-1),StartIndex_j:EndIndex_j) = Model.constraints{i}.Phi_qj(t, qi, []);

        index = index + length(Phi_c);
    elseif(Model.constraints{i}.body2 == 0) %Connected to ground
        StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
        EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
        qi = [q(StartIndex_i:EndIndex_i);[0;0;0;1;0;0;0]];

        Phi_c = Model.constraints{i}.Phi(t, qi, []);

        Phi(index:index+(length(Phi_c)-1)) = Phi_c;
        Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);

        index = index + length(Phi_c);                
    else
        StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
        EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
        StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
        EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
        qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];

        Phi_c = Model.constraints{i}.Phi(t, qi, []);

        Phi(index:index+(length(Phi_c)-1)) = Phi_c;
        Phi_q(index:index+(length(Phi_c)-1),StartIndex_i:EndIndex_i) = Model.constraints{i}.Phi_qi(t, qi, []);
        Phi_q(index:index+(length(Phi_c)-1),StartIndex_j:EndIndex_j) = Model.constraints{i}.Phi_qj(t, qi, []);

        index = index + length(Phi_c);
    end
end

assignin('base','Phi_q',Phi_q);


RHS = zeros(size(Phi_q,2),1);

for i = 1:length(Model.bodies)
    StartIndex_i = Model.bodies(i).StartIndex;
    
    RHS(StartIndex_i:StartIndex_i+2) = -Model.bodies(i).mass*eye(3) * Model.bodies(i).qdd(1:3,1);
    RHS(StartIndex_i+3:StartIndex_i+6) = -4*G(Model.bodies(i).q(4:7,1))'*diag(Model.bodies(i).MOI)*G(Model.bodies(i).q(4:7,1)) * Model.bodies(i).qdd(4:7,1) + ...
                                          8*G(Model.bodies(i).qd(4:7,1))'*diag(Model.bodies(i).MOI)*G(Model.bodies(i).qd(4:7,1)) * Model.bodies(i).q(4:7,1);
end

RHS

for i = 1:length(Model.forces)
    StartIndex_i = Model.bodies(Model.forces{i}.body1).StartIndex;
    
    RHS(StartIndex_i:StartIndex_i+6) = RHS(StartIndex_i:StartIndex_i+6) + Model.forces{i}.Q(t,Model.bodies(i).q(:,1),Model.bodies(i).qd(:,1));
end


RHS

lambda = Phi_q\RHS


end


function Model = read_Model(filename)
    %Model=loadjson(filename);
    Model = jsondecode(fileread(filename)); %Requires MATLAB R2016b
    
    %Save Body ID and Constrain ID vectors  
    Model.BodyIDs = zeros(length(Model.bodies),1);
    Model.ConstraintIDs = zeros(length(Model.constraints),1);
    Model.NumGeneralizedCoordinates = 7*ones(length(Model.bodies),1);
    
    %Process Defined Bodies
    for i = 1:length(Model.bodies)
        Model.bodies(i).q0 = Model.bodies(i).q0';
        Model.bodies(i).qd0 = Model.bodies(i).qd0';
        Model.BodyIDs(i) = Model.bodies(i).id;
    end
    
    Model.constraintDefinitions = Model.constraints;
    Model.constraints = {};
    
    
    %Process Defined Constraints
    for i = 1:length(Model.constraintDefinitions)
        
        %Change Body indices in constraints to match the index in which
        %they were defined.
        if(Model.constraintDefinitions(i).body1 ~= 0)
            Model.constraintDefinitions(i).body1 = find(Model.constraintDefinitions(i).body1 == Model.BodyIDs,1);
        end
        if(max(strcmpi('body2',fieldnames(Model.constraintDefinitions(i))))==1)
            if(Model.constraintDefinitions(i).body2 ~= 0)
                Model.constraintDefinitions(i).body2 = find(Model.constraintDefinitions(i).body2 == Model.BodyIDs,1);
            end
        end        
        
        switch Model.constraintDefinitions(i).type
            case 'Revolute'
                Model.constraints{i} = Constraint_Revolute(Model.constraintDefinitions(i));  
            otherwise
                error(['ERROR: Constraint type "',char(Model.constraintDefinitions(i).type),'" is not supported']);
        end
        Model.ConstraintIDs(i) = Model.constraintDefinitions(i).id;
    end   
    
    %Add in a Euler Parameter constraint for each body
    for i = 1:length(Model.bodies)
        Definition.name = ['P_',num2str(i)];
        Definition.id = length(Model.constraintDefinitions(i))+i;
        Definition.body1 = i;
        Model.constraints{end+1} = Constraint_P(Definition);  
        Model.ConstraintIDs(end+1) = Definition.id;
    end
    
    
    if(max(strcmpi(fieldnames(Model),'forces'))==0)
        Model.forces = {};
        Model.ForceIDs = [];
    end
    
    Model.forceDefinitions = Model.forces;
    Model.forces = {};
    
    %Process Defined Forces
    for i = 1:length(Model.forceDefinitions)
        
        %Change Body indices in forces to match the index in which
        %they were defined.
        Model.forceDefinitions(i).body1 = find(Model.forceDefinitions(i).body1 == Model.BodyIDs,1);
        
        switch Model.forces{i}.type
            case 'PointForce'
                Model.forces{i} = GenForce_PointForce(Model.forceDefinitions(i));
            otherwise
                error(['ERROR: Force type "',char(Model.forceDefinitions(i).type),'" is not supported']);
        end
        Model.ForceIDs(i) = Model.forces{i}.id;
    end
    NumUserDefinedForces = length(Model.forces);
    
    if(NumUserDefinedForces == 0)
        Model.forces = {};
    end
    
    %Add in Point forces for gravity
    for i = 1:length(Model.bodies)
        force.name = ['Gravity - Body',num2str(Model.BodyIDs(i))];
        force.id = max(Model.ForceIDs)+1;
        if(isempty(force.id))
            force.id = 1;
        end
        force.body1 = i;
        force.sP1 = [0;0;0];
        force.frame = 'GRF';
        force.funX = num2str(Model.gravity(1)*Model.bodies(i).mass);
        force.funY = num2str(Model.gravity(2)*Model.bodies(i).mass);
        force.funZ = num2str(Model.gravity(3)*Model.bodies(i).mass);
        
        Model.forces{NumUserDefinedForces + i} = GenForce_PointForce(force);
        Model.ForceIDs(NumUserDefinedForces + i) = force.id;
    end    
    
    %Save the start index for each body
    for i = 1:length(Model.bodies)
        if(i == 1)
            Model.bodies(i).StartIndex = 1;
        else
            Model.bodies(i).StartIndex = 1+sum(Model.NumGeneralizedCoordinates(1:i-1));
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