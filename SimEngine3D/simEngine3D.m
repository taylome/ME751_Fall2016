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
    model_name = [pathstr,'\models\me751_HW06.mdl'];
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
    Model = Run_KinematicAnalysis(Model,Settings);
elseif(strcmpi(Model.simulation.Type,'Dynamics'))
    error('Dynamic Analysis has not been written yet');
end
toc()

end

function Model = Run_KinematicAnalysis(Model,Settings)

NumTimeSteps = Model.simulation.tend/Model.simulation.stepSize+1;

for i = 1:length(Model.bodies)
     Model.bodies(i).q = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
     Model.bodies(i).qd = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
     Model.bodies(i).qdd = zeros(Model.NumGeneralizedCoordinates(i),NumTimeSteps);
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
%        guess - Return to Step 2 until final time is reached
%--------------------------------------------------------------------------
    for i = 1:length(Model.bodies)
        Model.bodies(i).q(:,tindex) = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        Model.bodies(i).qd(:,tindex) = qd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        Model.bodies(i).qdd(:,tindex) = qdd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
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
    
    assignin('base','Model',Model);
    
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
    
    
    %Save the start index for each body
    for i = 1:length(Model.bodies)
        if(i == 1)
            Model.bodies(i).StartIndex = 1;
        else
            Model.bodies(i).StartIndex = 1+sum(Model.NumGeneralizedCoordinates(1:i-1));
        end
    end    
end