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
    model_name = [pathstr,'\models\me751_HW08P2.mdl'];
end

%--------------------------------------------------------------------------
% Settings
%--------------------------------------------------------------------------
Settings.NewtonRaphsonMaxIterations = 50;
%Settings.NewtonRaphsonTol = 1E-8;
Settings.NewtonRaphsonTol = 1E-4;
Settings.model_name = model_name;
%--------------------------------------------------------------------------
 
disp(['Reading Model: ',model_name]);

[Model] = read_Model(model_name);
% assignin('base','Model',Model)


if(strcmpi(Model.simulation.Type,'Kinematics'))
    disp('Executing Kinematic Analysis...');
    Model = Run_KinematicAnalysis(Model,Settings,false);
elseif(strcmpi(Model.simulation.Type,'Inverse Dynamics'))
    disp('Executing Inverse Dynamics Analysis...');
    Model = Run_KinematicAnalysis(Model,Settings,true);    
elseif(strcmpi(Model.simulation.Type,'Dynamics'))
    disp('Executing Dynamic Analysis...');
    Model = Run_DynamicsAnalysis(Model,Settings);   
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
        [Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q);
        
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
    Nu = calc_Nu(Model,t,q);
    qd = Phi_q\Nu;

%--------------------------------------------------------------------------
%Step 4: Acceleration Analysis - Solve Phi_q*qdot = Gamma
%--------------------------------------------------------------------------
    Gamma = calc_Gamma(Model,t,q,qd);    
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
        lambda = Phi_q'\RHS;
        
        Model.lambda(:,tindex) = lambda;
        
        % Solve for the reaction forces and torques for each constraint
        Model = Save_ReactionForceTorque(Model,q,Phi_q,lambda,tindex);
    end
    
end

Model.time = ((0:NumTimeSteps-1)*Model.simulation.stepSize)';


end

function Model = Run_DynamicsAnalysis(Model,Settings)

NumTimeSteps = Model.simulation.tend/Model.simulation.stepSize+1;
OutputTimePoints = linspace(0,Model.simulation.tend,Model.simulation.outputSteps+1);
OutputTimePointsCounter = 2;


for i = 1:length(Model.bodies)
     Model.bodies(i).q = zeros(Model.NumGeneralizedCoordinates(i),Model.simulation.outputSteps);
     Model.bodies(i).qd = zeros(Model.NumGeneralizedCoordinates(i),Model.simulation.outputSteps);
     Model.bodies(i).qdd = zeros(Model.NumGeneralizedCoordinates(i),Model.simulation.outputSteps);
end

for i = 1:length(Model.constraints)
    Model.ConstraintReactions(i).Body1 = zeros(6,Model.simulation.outputSteps);
    Model.ConstraintReactions(i).Body2 = zeros(6,Model.simulation.outputSteps);
end

%--------------------------------------------------------------------------
% Solve for the generalized accelerations and Lagrange multipliers at t = 0
%--------------------------------------------------------------------------
t = 0;

q = zeros(sum(Model.NumGeneralizedCoordinates),1);
qd = zeros(sum(Model.NumGeneralizedCoordinates),1);
for i = 1:length(Model.bodies)
    q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i)) = Model.bodies(i).q0;
    qd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i)) = Model.bodies(i).qd0;
end


%Solve for the LHS Matrix
MassJPMatrix = calc_MassJPMatrix(Model,q);

[Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q);

        
%Solve for the RHS Vector
[QA] = calc_QA(Model,t,q,qd);

[Gamma] = calc_Gamma(Model,t,q,qd);

%Now solve for qdd and lambda
qdd_lambda = [MassJPMatrix,Phi_q';Phi_q,zeros(size(Phi_q,1))]\[QA;Gamma];
qdd = qdd_lambda(1:sum(Model.NumGeneralizedCoordinates));
lambda = qdd_lambda(sum(Model.NumGeneralizedCoordinates)+1:end);

%Now Save the results for t=0 back to the model structure
for i = 1:length(Model.bodies)
    Model.bodies(i).q(:,1) = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
    Model.bodies(i).qd(:,1) = qd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
    Model.bodies(i).qdd(:,1) = qdd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
end

Model.lambda = zeros(length(lambda),Model.simulation.outputSteps);
Model.lambda(:,1) = lambda;
Model.NR = zeros(1,Model.simulation.outputSteps);
Model.correction_norm = zeros(1,Model.simulation.outputSteps);
Model.correction = zeros(length(qdd_lambda),Model.simulation.outputSteps);
Model.CondPsi = zeros(1,Model.simulation.outputSteps);

%Now Solve for the reaction forces/torques for each of the defined
%constraints
Model = Save_ReactionForceTorque(Model,q,Phi_q,lambda,1);


% assignin('base','Model',Model);
% return

%Check to make sure we have a healthy starting point
Beta0 = 1;
h = Model.simulation.stepSize;
MassJPMatrix = calc_MassJPMatrix(Model,q);

[Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q);

Psi = [MassJPMatrix,Phi_q';Phi_q,zeros(size(Phi_q,1))];       

%Solve for the Applied Force/Torque Vector
[QA] = calc_QA(Model,t,q,qd);       

%Now solve for the correction to qdd_lambda
g = [MassJPMatrix*qdd+Phi_q'*lambda-QA;...
    1/(Beta0^2*h^2)*Phi];

if(max(g)>Settings.NewtonRaphsonTol)
    msgbox('Check Model Definition File.  Initial solution may not be valid');
end

%--------------------------------------------------------------------------
% Now that we've got a good starting point, start into the simulation loop
%--------------------------------------------------------------------------

%Allocate Placeholder for up to BDF of order 6
q_prev = repmat(q,1,6);
qd_prev = repmat(qd,1,6);
qdd_prev = repmat(qdd,1,6);
lambda_prev = repmat(lambda,1,6);

% qdd = qdd*0;
% lambda = lambda*0;
%Model.simulation.Solver = 'BDF1';
h = Model.simulation.stepSize;
for tindex = 2:NumTimeSteps
    
    %Shift the values from the last iteration onto the previous values
    %stacks
    q_prev = [q,q_prev(:,1:end-1)];
    qd_prev = [qd,qd_prev(:,1:end-1)];
    qdd_prev = [qdd,qdd_prev(:,1:end-1)];
    lambda_prev = [lambda,lambda_prev(:,1:end-1)];
    
    t = Model.simulation.stepSize * (tindex-1);
    
    %Start the Quasi-NR iterations.    
    for NR = 1:Settings.NewtonRaphsonMaxIterations
        if((tindex<3)||(strcmpi(Model.simulation.Solver,'BDF1')))
            qd = qd_prev(:,1)+h*qdd;
            q = q_prev(:,1)+h*qd;
            Beta0 = 1;
            BDForder = 1;
        else %Assume BDF order 2 for now, add other orders later
            qd = 4/3*qd_prev(:,1)-1/3*qd_prev(:,2)+2/3*h*qdd;
            q = 4/3*q_prev(:,1)-1/3*q_prev(:,2)+2/3*h*qd;
            Beta0 = 2/3;
            BDForder = 2;
        end
                    
        MassJPMatrix = calc_MassJPMatrix(Model,q);

        [Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q);
            
        %Solve for the Applied Force/Torque Vector
        [QA] = calc_QA(Model,t,q,qd);       
        
        %Now solve for the correction to qdd_lambda
        g = [MassJPMatrix*qdd+Phi_q'*lambda-QA;...
            1/(Beta0^2*h^2)*Phi];
        
        switch(lower(Model.simulation.NRMethod))
            case {lower('QuasiFullNR'),lower('QuasiModifiedNR')}
                if((NR == 1) || strcmpi(Model.simulation.NRMethod,'QuasiFullNR'))
                    Psi = [MassJPMatrix,Phi_q';Phi_q,zeros(size(Phi_q,1))];
                end
            case {lower('NumericFullNR'),lower('NumericModifiedNR')}
                if((NR == 1)||strcmpi(Model.simulation.NRMethod,'NumericFullNR'))
                    Psi = zeros(length(g));
                    delta = 0.0001;
                    for i = 1:length(qdd_lambda)
                        qdd_lambda_delta = 0*qdd_lambda;
                        qdd_lambda_delta(i) = delta;
                        qdd_lambda_temp = qdd_lambda + qdd_lambda_delta;
                        qdd_temp = qdd_lambda_temp(1:sum(Model.NumGeneralizedCoordinates));
                        lambda_temp = qdd_lambda_temp(sum(Model.NumGeneralizedCoordinates)+1:end);
                        Psi(:,i) = (calc_g(Model,t,q_prev,qd_prev,qdd_temp,lambda_temp,h,BDForder)-g)/delta;
                    end
                end
            case {lower('FullNR'),lower('ModifiedNR')}
                if((NR == 1)||strcmpi(Model.simulation.NRMethod,'FullNR'))
                    Psi = [MassJPMatrix,Phi_q';Phi_q,zeros(size(Phi_q,1))];

                    %Add in h^2*Bo^2*[Jp(pn)pdd]_p
                    for i = 1:length(Model.bodies)
                        qi = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
                        qddi = qdd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
                        J = diag(Model.bodies(i).MOI);

                        Psi((4+7*(i-1)):(7+7*(i-1)),(4+7*(i-1)):(7+7*(i-1))) = Psi((4+7*(i-1)):(7+7*(i-1)),(4+7*(i-1)):(7+7*(i-1))) + ...
                            (Beta0^2*h^2)*Jp_p_pdd_p(J,qi(4:7),qddi(4:7));
                    end

                    %Add in Phi_q_q terms
                    index = 1;
                    for i = 1:length(Model.constraints)
                        if(Model.constraints{i}.numbodies == 1)
                            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                            qi = q(StartIndex_i:EndIndex_i);

                            Phi_c = Model.constraints{i}.Phi(t, qi, []);

                            lambdai = lambda(index:index-1+length(Phi_c));

                            Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) = Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qi_lambda_qi(qi,lambdai);

                            index = index + length(Phi_c);
                        elseif(Model.constraints{i}.body1 == 0) %Connected to ground
                            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                            qi = [[0;0;0;1;0;0;0];q(StartIndex_j:EndIndex_j)];

                            Phi_c = Model.constraints{i}.Phi(t, qi, []);

                            lambdai = lambda(index:index-1+length(Phi_c));

                            Psi(StartIndex_j:EndIndex_j,StartIndex_j:EndIndex_j) = Psi(StartIndex_j:EndIndex_j,StartIndex_j:EndIndex_j) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qj_lambda_qj(qi,lambdai);

                            index = index + length(Phi_c);
                        elseif(Model.constraints{i}.body2 == 0) %Connected to ground
                            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                            qi = [q(StartIndex_i:EndIndex_i);[0;0;0;1;0;0;0]];

                            Phi_c = Model.constraints{i}.Phi(t, qi, []);

                            lambdai = lambda(index:index-1+length(Phi_c));

                            Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) = Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qi_lambda_qi(qi,lambdai);

                            index = index + length(Phi_c);                
                        else
                            StartIndex_i = Model.bodies(Model.constraints{i}.body1).StartIndex;
                            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body1);
                            StartIndex_j = Model.bodies(Model.constraints{i}.body2).StartIndex;
                            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.constraints{i}.body2);
                            qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];

                            Phi_c = Model.constraints{i}.Phi(t, qi, []);

                            lambdai = lambda(index:index-1+length(Phi_c));

                            Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) = Psi(StartIndex_i:EndIndex_i,StartIndex_i:EndIndex_i) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qi_lambda_qi(qi,lambdai);
                            Psi(StartIndex_i:EndIndex_i,StartIndex_j:EndIndex_j) = Psi(StartIndex_i:EndIndex_i,StartIndex_j:EndIndex_j) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qi_lambda_qj(qi,lambdai);
                            Psi(StartIndex_j:EndIndex_j,StartIndex_i:EndIndex_i) = Psi(StartIndex_j:EndIndex_j,StartIndex_i:EndIndex_i) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qj_lambda_qi(qi,lambdai);
                            Psi(StartIndex_j:EndIndex_j,StartIndex_j:EndIndex_j) = Psi(StartIndex_j:EndIndex_j,StartIndex_j:EndIndex_j) + (Beta0^2*h^2)*Model.constraints{i}.Phi_qj_lambda_qj(qi,lambdai);

                            index = index + length(Phi_c);
                        end
                    end

                    %Add in Force Terms (numerically calculate derivatives
                    %since functions might be unknown.
                    QA_q = zeros(size(QA));
                    QA_qd = zeros(size(QA));
                    delta = 0.0001;
                    for i = 1:length(q)
                        q_delta = 0*q;
                        q_delta(i) = delta;
                        QA_q(:,i) = (calc_QA(Model,t,q+q_delta,qd)-QA)/delta;
                        QA_qd(:,i) = (calc_QA(Model,t,q,qd+q_delta)-QA)/delta;
                    end
                    Psi(1:length(q),1:length(q)) = Psi(1:length(q),1:length(q)) - (Beta0^2*h^2)*(QA_q+QA_qd);
                end
            otherwise
                error('Unknown NR Method Type');
        end
                
            
        correction = Psi\-g;        
        
        %Now solve for qdd and lambda
        qdd_lambda = qdd_lambda+correction;
        qdd = qdd_lambda(1:sum(Model.NumGeneralizedCoordinates));
        lambda = qdd_lambda(sum(Model.NumGeneralizedCoordinates)+1:end);
    
        %stop = norm(correction);
        %if(max(abs(correction))<Settings.NewtonRaphsonTol)
        if(norm(correction)<Settings.NewtonRaphsonTol)
            break
        end
    end
    
    %Calcuate the final update to q & qd based on qdd
    if((tindex<3)||(strcmpi(Model.simulation.Solver,'BDF1')))
        qd = qd_prev(:,1)+h*qdd;
        q = q_prev(:,1)+h*qd;
        Beta0 = 1;
    else %Assume BDF order 2 for now, add other orders later
        qd = 4/3*qd_prev(:,1)-1/3*qd_prev(:,2)+2/3*h*qdd;
        q = 4/3*q_prev(:,1)-1/3*q_prev(:,2)+2/3*h*qd;
        Beta0 = 2/3;
    end
        
%--------------------------------------------------------------------------
% Store the q vectors into the bodies if a print point is reached
%--------------------------------------------------------------------------
    if(t>=OutputTimePoints(OutputTimePointsCounter))
        for i = 1:length(Model.bodies)
            Model.bodies(i).q(:,OutputTimePointsCounter) = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
            Model.bodies(i).qd(:,OutputTimePointsCounter) = qd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
            Model.bodies(i).qdd(:,OutputTimePointsCounter) = qdd(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        end
        

        %Now Solve for the reaction forces/torques for each of the defined
        %constraints
        Model = Save_ReactionForceTorque(Model,q,Phi_q,lambda,OutputTimePointsCounter);

        Model.lambda(:,OutputTimePointsCounter) = lambda;
        Model.NR(:,OutputTimePointsCounter) = NR;
        Model.Phi(:,OutputTimePointsCounter) = Phi;
        Model.PhiMax(:,OutputTimePointsCounter) = max(Phi);
        Model.g(:,OutputTimePointsCounter) = g;
        Model.gMax(:,OutputTimePointsCounter) = max(g);
        Model.QA(:,OutputTimePointsCounter) = QA;
        Model.correction_norm(OutputTimePointsCounter) = norm(correction);
        Model.correction(:,OutputTimePointsCounter) = correction;
        Model.CondPsi(OutputTimePointsCounter) = cond(Psi);
        
        OutputTimePointsCounter = OutputTimePointsCounter + 1;
        if(mod(OutputTimePointsCounter,100)==1)
            fprintf(['t = ',num2str(t),' s - ',num2str(tindex/NumTimeSteps*100),'%% complete\n']);
        end
    end
    
end
Model.time = ((0:NumTimeSteps-1)*Model.simulation.stepSize)';

assignin('base','Model',Model)


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
            case 'Spherical'
                Model.constraints{i} = Constraint_Spherical(Model.constraintDefinitions(i)); 
            otherwise
                error(['ERROR: Constraint type "',char(Model.constraintDefinitions(i).type),'" is not supported']);
        end
        Model.ConstraintIDs(i) = Model.constraintDefinitions(i).id;
    end   
    
    %Add in a Euler Parameter constraint for each body
    for i = 1:length(Model.bodies)
        Definition.name = ['P_',num2str(i)];
        Definition.id = length(Model.constraintDefinitions)+i;
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


function MassJPMatrix = calc_MassJPMatrix(Model,q)
    MassJPMatrix = eye(length(Model.bodies)*7);
    for i = 1:length(Model.bodies)
        qi = q(Model.bodies(i).StartIndex:Model.bodies(i).StartIndex-1+Model.NumGeneralizedCoordinates(i));
        
        MassJPMatrix((1+7*(i-1)):(3+7*(i-1)),(1+7*(i-1)):(3+7*(i-1))) = Model.bodies(i).mass*eye(3);
        MassJPMatrix((4+7*(i-1)):(7+7*(i-1)),(4+7*(i-1)):(7+7*(i-1))) = 4*G(qi(4:7))'*diag(Model.bodies(i).MOI)*G(qi(4:7));
    end
end

function [Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q)
    Phi_q = zeros(size(q,1));
    Phi = zeros(size(q));
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
    Phi_q(index:end,:)=[]; %Eliminate any missing constraint equations
    Phi(index:end)=[]; %Eliminate any missing constraint equations
end

function Nu = calc_Nu(Model,t,q)
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
    Nu(index:end,:)=[]; %Eliminate any missing constraint equations
end

function Gamma = calc_Gamma(Model,t,q,qd)
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
    Gamma(index:end,:)=[]; %Eliminate any missing constraint equations
end

function QA = calc_QA(Model,t,q,qd)
    QA = zeros(size(q));
    for i = 1:length(Model.forces)
        if(Model.forces{i}.numbodies == 1)
            StartIndex_i = Model.bodies(Model.forces{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.forces{i}.body1);
            qi = q(StartIndex_i:EndIndex_i);
            qdi = qd(StartIndex_i:EndIndex_i);

            QA(StartIndex_i:EndIndex_i) = QA(StartIndex_i:EndIndex_i) + Model.forces{i}.Q(t, qi, qdi);
        else
            StartIndex_i = Model.bodies(Model.forces{i}.body1).StartIndex;
            EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.forces{i}.body1);
            StartIndex_j = Model.bodies(Model.forces{i}.body2).StartIndex;
            EndIndex_j = StartIndex_j-1+Model.NumGeneralizedCoordinates(Model.forces{i}.body2);
            qi = [q(StartIndex_i:EndIndex_i);q(StartIndex_j:EndIndex_j)];
            qdi = [qd(StartIndex_i:EndIndex_i);qd(StartIndex_j:EndIndex_j)];

            QA(StartIndex_i:EndIndex_i) = QA(StartIndex_i:EndIndex_i) + Model.forces{i}.Q1(t, qi, qdi);
            QA(StartIndex_j:EndIndex_j) = QA(StartIndex_i:EndIndex_i) + Model.forces{i}.Q2(t, qi, qdi);
        end
    end

    for i = 1:length(Model.bodies)
        StartIndex_i = Model.bodies(Model.forces{i}.body1).StartIndex;
        EndIndex_i = StartIndex_i-1+Model.NumGeneralizedCoordinates(Model.forces{i}.body1);
        qi = q(StartIndex_i:EndIndex_i);
        qdi = qd(StartIndex_i:EndIndex_i);

        QA(StartIndex_i+3:EndIndex_i) = QA(StartIndex_i+3:EndIndex_i) + 8*G(qdi(4:7))'*diag(Model.bodies(i).MOI)*G(qdi(4:7)) * qi(4:7);
    end
end

function g = calc_g(Model,t,q_prev,qd_prev,qdd,lambda,h,BDForder)
        if(BDForder == 1)
            qd = qd_prev(:,1)+h*qdd;
            q = q_prev(:,1)+h*qd;
            Beta0 = 1;
        else %Assume BDF order 2 for now, add other orders later
            qd = 4/3*qd_prev(:,1)-1/3*qd_prev(:,2)+2/3*h*qdd;
            q = 4/3*q_prev(:,1)-1/3*q_prev(:,2)+2/3*h*qd;
            Beta0 = 2/3;
        end
        
        MassJPMatrix = calc_MassJPMatrix(Model,q);

        [Phi,Phi_q] = calc_Phi_Phi_q(Model,t,q);      

        %Solve for the Applied Force/Torque Vector
        [QA] = calc_QA(Model,t,q,qd);       
        
        %Now solve for the correction to qdd_lambda
        g = [MassJPMatrix*qdd+Phi_q'*lambda-QA;...
            1/(Beta0^2*h^2)*Phi];
end

function Model = Save_ReactionForceTorque(Model,q,Phi_q,lambda,Save_Index)
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
            currentlambda = lambda(index + c-1);
            Force = -Phi_qj(c,1:3)'*currentlambda; %Global Frame
            Torque = -1/2*G(qj(4:7))*Phi_qj(c,4:7)'*currentlambda-Tilde(Model.constraints{i}.marker2.origin)*A(qj(4:7))'*Force;
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
            currentlambda = lambda(index + c-1);
            Force = -Phi_qi(c,1:3)'*currentlambda; %Global Frame
            Torque = -1/2*G(qi(4:7))*Phi_qi(c,4:7)'*currentlambda-Tilde(Model.constraints{i}.marker1.origin)*A(qi(4:7))'*Force;
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
            currentlambda = lambda(index + c-1);
            Force = -Phi_qi(c,1:3)'*currentlambda; %Global Frame
            Torque = -1/2*G(qi(4:7))*Phi_qi(c,4:7)'*currentlambda-Tilde(Model.constraints{i}.marker1.origin)*A(qi(4:7))'*Force;
            Torque = A(qi(4:7))*Torque; %Local -> Global Frame

            AccumForceTrq_Bodyi = AccumForceTrq_Bodyi + [Force;Torque];

            Force = -Phi_qj(c,1:3)'*currentlambda; %Global Frame
            Torque = -1/2*G(qj(4:7))*Phi_qj(c,4:7)'*currentlambda-Tilde(Model.constraints{i}.marker2.origin)*A(qj(4:7))'*Force;
            Torque = A(qj(4:7))*Torque; %Local -> Global Frame

            AccumForceTrq_Bodyj = AccumForceTrq_Bodyj + [Force;Torque];                    
        end                

    end            

    Model.ConstraintReactions(i).Body1(:,Save_Index) = AccumForceTrq_Bodyi;
    Model.ConstraintReactions(i).Body2(:,Save_Index) = AccumForceTrq_Bodyj;

    index = index + NumEquations;
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
function mess = Jp_p_pdd_p(J,p,pdd)
    a = J*G(p)*pdd;
    mess = 4*(T(a)-G(p)'*J*G(pdd));
end
