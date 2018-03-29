classdef AnkleModel < handle
    % Simple model of standing postural stability, consisting of foot and
    % body segments, and two muscles that create moments about the ankles,
    % tibialis anterior. 
    
    methods (Static)
        
        function result = tibialisLength(theta)
            % Soleus length as a function of ankle angle. 
            % 
            % theta: body angle (up from prone horizontal)
            % result: tibialis anterior length
            
            origin = zeros(2,length(theta));
            for i = 1:length(theta)
                origin(:,i) = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]*[.3; -.03];  
            end
            insertion = [.06; -.03];
            difference = origin - insertion;
            result = (difference(1,:).^2 + difference(2,:).^2).^.5;

            % return result with same shape as theta
            if size(theta, 1) > size(theta, 2)
                result = result';
            end
        end        
        
        function simulate(control, T)
            % Runs a simulation of the model and plots results. 
            % 
            % control: 0 means no control law should be usd to stabilize
            %   the model; 1 means a control law should be used
            % T: total time to simulate, in seconds
            
            % TO DO: set parameters
            % plot AE
        end        
    end
end

function dx_dt = dynamics(t, x, TA, control) 
    % Right-hand side of the dynamic equation of the model. 
    % 
    % t: time of current step (s)
    % x: model state: [ankle angle, angular velocity, soleus CE length, TA
    %   CE length]
    % S: soleus muscle object (a HillTypeModel)
    % TA: tibialis anterior muscle object (a HillTypeModel)
    
    % TO DO: define original params

    taLengthSE = TA.getNormalizedLengthSE(StabilityModel.tibialisLength(x(1)),x(4));
    
    % mTA: moment created by TA
    mTA = TA.f0M * TA.forceLengthSE( taLengthSE )* dTA;
    
    % a_dt: change in angle over time
    % v_dt: change in angular velocity over time
    % ta_dt: change in TA CE length over time
    
    % TO DO: finish dynamics eqn.
    a_dt = x(2);
    v_dt = 
    ta_dt =  TA.getVelocity(aTA, x(4), taLengthSE);
  
    dx_dt = [a_dt v_dt ta_dt]';
end

function result = getGravityMoment(angle)
    % angle: angle of body segment (up from horizontal)
    % result: moment about ankle due to force of gravity on body
    
%     m = 75; % body segment mass 
%     lCOM = 1; % distance from ankle to body segment centre of mass
%     g = 9.81; 
%     result = m*g*lCOM*sin(angle-pi/2);

    % TO DO: calculate gravity moment of ankle
end

function result = getActivationTA(angle)
    % TO DO: get activation of TA at any given point in swing phase
end

function result = getMomentArm(angle)
    % TO DO: look at angle of shank

end
