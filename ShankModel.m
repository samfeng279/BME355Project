classdef ShankModel < handle
    % Simple model of the shank, represented by a pendulum
    
    methods (Static)
        % TO DO: pendulum code
        
        function simulate(control, T)
            % TO DO: simulate and plot graphs
        end        
    end
end

function dx_dt = dynamics() 
    % TO DO: model pendulum at ankle joint at any given time
end

function result = getGravityMoment(angle)
    % angle: angle of body segment (up from horizontal)
    % result: moment about ankle due to force of gravity on body
    
    % TO DO: get gravity moment about shank (and ankle?)
end

function result = getShankAngle()
    % TO DO: get shank angle, will affect ankle moment arm 
end

function result = getInertialForces()
    % TO DO: get inertial forces as due to swinging of ankle
end