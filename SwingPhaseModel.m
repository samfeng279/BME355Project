classdef SwingPhaseModel < handle
    % Simple model of standing postural stability, consisting of foot and
    % body segments, and two muscles that create moments about the ankles,
    % tibialis anterior and soleus. 
    
    properties (Access = public)
        shankAngle = 5.926; % Initial angle of shank from -pi/2
        footAngle = 6.574; % Initial angle of COM of foot from -pi/2
        shankAngVel = 2.042;  % Initial angular velocity of shank, stimulated
        footAngVel = 2.460; % Initial angular velocity of COM of foot, stimulated
        shankAngVelDrop = 2.042*0.47; % Initial angular velocity of shank, unstimulated
        footAngVelDrop = 2.460*0.47; % Initial angular velocity of COM of foot, unstimulated
        shankMass = 0.0433*80.7; 
        footMass = 0.0137*80.7;  
        shankLength = 0.4318; 
        ankleFootLength = 0.07136; % Distance from ankle to COM of foot
        COMToeLength = 0.117883; % Distance from COM of foot to toes
        duration = 0.65; % Duration of simulation
        fps = 100; % Frames per sec - changes resolution of animation 
        modelledTA; % TA activation based off regression, array of TA values and time
        pointNum = 5000; % Number of points in modelledTA
        gaitInterval = 1.41; % Time required for one gait cycle
        swingInterval = [0.6 1]; % percentage of the gait cycle that is the swing phase
        frequency = 34; % Ideal frequency of stimulation
    end
    
    methods (Static)
        
        function result = soleusLength(theta)
            % Soleus length as a function of ankle angle. 
            % 
            % theta: body angle (up from prone horizontal)
            % result: soleus length
            
            origin = zeros(2,length(theta));
            for i = 1:length(theta) %this loop optionally handles a list of theta 
                origin(:,i) = [cos(theta(i)) -sin(theta(i)); sin(theta(i)) cos(theta(i))]*[.3; .03];                
            end
            insertion = [-.05; -.02]; 
            difference = origin - insertion;
            result = (difference(1,:).^2 + difference(2,:).^2).^.5;
            
            % return result with same shape as theta
            if size(theta, 1) > size(theta, 2)
                result = result';
            end
        end
        
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
        
        % function for running regression on the muscle activation graph to
        % get modellig equations
        function result = getActivationRegression(data, gaitInterval)
            % updating time values to fit with the gait interval
            for i = 1:size(data,1)
                data(i,1) = data(i,1)*(10^-2)*gaitInterval;
            end

            % finding the regressed line for each point to map the curve
            delta = 5;
            regressions = []; % matrix of time points and their corresponding coeffs
            for i = 1:size(data,1)
                % finding start and end points of data interval being regressed
                if i < delta + 1
                    s = 1;
                    e = i + delta;
                elseif (size(data,1) - i) < delta
                    s = i - delta;
                    e = size(data,1);
                else
                    s = i - delta;
                    e = i + delta;
                end

                %take care of a zero case, otherwise we only begin when t > 0 
                if i == 1
                    x = [0; data(1:9,1)];
                    y = [data(1,2); data(1:9,2)];
                    coeffs = polyfit(x,y,2);
                    regressions = [regressions; 0 coeffs(1) coeffs(2) coeffs(3)];
                end

                x = data(s:e,1);
                y = data(s:e,2);
                coeffs = polyfit(x,y,2);
                regressions = [regressions; data(i,1) coeffs(1) coeffs(2) coeffs(3)];
            end
            
            result = regressions;            
        end
        
        % fuction for obtaining simulated activation values from modelled
        % curve
        function result = getSimulatedActivations(regressions, gaitInterval, swingInterval)
            % data set gives points to model the activation of TA, and from
            % here we calculate activations of TA through the gait cycle
            % using the regression -> create new time array (select point #)
            times = linspace(0, gaitInterval, 5000);
            simulated = [];
            for i = 1:size(times,2)
                for j = 1:size(regressions,1)-1
                    % finding which interval of the regression data the new time falls into
                    if (times(i) >= regressions(j,1)) && (times(i) <= regressions(j+1,1))
                        % now that we've found the interval, we can approx the
                        % activation value using either j's or j+1's polynomial
                        a = polyval([regressions(j, 2) regressions(j, 3) regressions(j, 4)], times(i));
                        simulated = [simulated; times(i) a];
                        break;
                    end
                end
            end
            
            % clipping data for the swing phase using swingInterval
            minIndex = round(swingInterval(1)*size(simulated,1));
            maxIndex = round(swingInterval(2)*size(simulated,1));
            simulated = simulated(minIndex:maxIndex,1:2);
            result = simulated;
        end
        
    end
    
    methods (Access = public)
        
        %20-50Hz
        % will be loading data = csvread('TA_STIM.csv'); 
        function getModelledTA(SP, data)
            regressions = SwingPhaseModel.getActivationRegression(data, SP.gaitInterval);
            simulatedTA = SwingPhaseModel.getSimulatedActivations(regressions, SP.gaitInterval, SP.swingInterval);
            
            % using frequency to find how many points exist between zeros and peaks
            delta = round((SP.pointNum/SP.gaitInterval)/SP.frequency);
            % identifying where peaks and zeros occur, as well as identifying the linear
            % relation between adjacents
            latestPeak = 1; latestZero = 1;
            peaks = []; zeros = [1];
            eqns = []; % store info in format: start_index end_index coeff(1) coeff(2)
            for i = 1:size(simulatedTA,1)    
                % identify where and when peaks should occur (to create jagged activation)
                if (latestPeak == 1) && ((i - latestPeak) > delta/2)
                    latestPeak = i;
                    peaks = [peaks; i];
                else
                    if (i - latestPeak) > delta
                        latestPeak = i;
                        peaks = [peaks; i];
                    elseif (i - latestZero) > delta
                        latestZero = i;
                        zeros = [zeros; i];
                    end  
                end

                % looking at the ascending curve to peak
                if latestZero < latestPeak
                    % check for empty eqns (check is not necessary for descending case
                    % as it is assumed ascending comes first
                    if isempty(eqns)
                        coeffs = polyfit([simulatedTA(latestZero,1) simulatedTA(latestPeak,1)], [simulatedTA(latestZero,2) simulatedTA(latestPeak,2)], 1);
                        eqns = [eqns; latestZero latestPeak coeffs(1) coeffs(2)];
                    % checking if the coeffs have already been found for these
                    % particular points
                    elseif ~any(eqns(:,1) == latestZero)
                        coeffs = polyfit([simulatedTA(latestZero,1) simulatedTA(latestPeak,1)], [0 simulatedTA(latestPeak,2)], 1);
                        eqns = [eqns; latestZero latestPeak coeffs(1) coeffs(2)];
                    end
                % looking at the descending curve from peak 
                elseif latestZero > latestPeak
                    if ~any(eqns(:,1) == latestPeak)
                        coeffs = polyfit([simulatedTA(latestPeak,1) simulatedTA(latestZero,1)], [simulatedTA(latestPeak,2) 0], 1);
                        eqns = [eqns; latestPeak latestZero coeffs(1) coeffs(2)];
                    end
                end
            end

            % finded modelled activation values
            activations = [];
            for i = 1:size(simulatedTA,1)
                if any(zeros(:) == i)
                    % if labelled as a zero, input into the activation matrix with
                    % activation of 0
                    activations = [activations; simulatedTA(i,1) 0];
                elseif any(peaks(:) == i)
                    % if labelled as a peak, input into the activation matrix with
                    % sample activation 
                    activations = [activations; simulatedTA(i,1) simulatedTA(i,2)];
                else
                    % else, based on the location of the point, use the coefficients found 
                    % to find the corresponding activation
                    % find new value based on linearity between adjacent zero and peak
                    for j = 1:size(eqns,1)
                        if  i < eqns(j,2) && i > eqns(j,1)
                            a = polyval([eqns(j,3) eqns(j,4)],simulatedTA(i,1));
                            activations = [activations; simulatedTA(i,1) a];
                        end
                    end
                end
            end

            SP.modelledTA = activations;
        end   
                
        function aTA = getActivationTA(SP, t)
            % Get activation of TA at a given time of step cycle
            
            % Find closest value of SP.modelledTA time to t 
            % First time of modelledTA considered is 0.8433 (offset)            
            [c, index] = min(abs(SP.modelledTA(:,1)'-t-0.8433));
            
            % Activation is the activation at that index
            aTA = SP.modelledTA(index,2);
        end
        
        function aS = getActivationS(SP, t)
           	% Activation of soleus constantly 0 for this model
            aS = 0;
        end
        
        function simulate(SP, stimulation)
            % If stimulation == 1, external activation applied to TA 
            % Else, no stimulation applied
            
            % Read in TA activation data 
            taData = csvread('TA_STIM.csv');
            % Regress data, returning array of TA values and time
            SP.getModelledTA(taData); 
            
            clc; figure;
            
            % Initialize TA and soleus muscles
            restLengthTA = SP.tibialisLength(pi/2);
            TA = HillTypeMuscle(2000, .6*restLengthTA, .4*restLengthTA);
            restLengthS = SP.soleusLength(pi/2);
            S = HillTypeMuscle(16000, .6*restLengthS, .4*restLengthS);
           
            % Angular velocities differ if foot is stimulated vs not
            if stimulation
                footVel = SP.footAngVel;
                shankVel = SP.shankAngVel;
            else                
                footVel = SP.footAngVelDrop;
                shankVel = SP.shankAngVelDrop;
            end
            
            % Set initial parameters for ODE
            ivp=[SP.shankAngle; SP.footAngle; shankVel; footVel; 1.045; 0.8129];
            
            % Simulate double pendulum behaviour
            SP.double_pendulum(ivp, SP.duration, SP.fps, TA, S, stimulation);
            
            %Simulate isometric forces
            SP.isometric(TA, 2000,.6*restLengthTA,.4*restLengthTA, stimulation);
        end 
        
        function xdot = dynamics(SP, t, x, TA, S, stimulation)
            % x: [shankAngle; footAngle; shankAngVel; footAngVel; CE TA length; CE soleus length];
            
            xdot = zeros(6,1);
            m1 = SP.shankMass; 
            m2 = SP.footMass; 
            l1 = SP.shankLength; 
            d1 = 0.494 * l1;
            l2 = SP.ankleFootLength;
            g = 9.81;
            I1 = 90;
            I2 = 80;
            I2COM = 85;
            
            % TA activation different if TA is stimulated  
            if stimulation
                aTA = 10000*SP.getActivationTA(t);
            else
                aTA = 0;
            end
            
            aS = SP.getActivationS(t);
            
            % Calculate moment about ankle due to TA
            taLength = SwingPhaseModel.tibialisLength(pi-x(1)-x(2));
            dTA = .03;
            fTA = 2000;
            taMoment = fTA*TA.forceLengthSE(TA.getNormalizedLengthSE(taLength,x(5)))*dTA;
            
            % Calculate moment about ankle due to soleus
            sLength = SwingPhaseModel.soleusLength(pi-x(1)-x(2));
            dS = .05;
            fS = 16000;
            sMoment = fS*S.forceLengthSE(S.getNormalizedLengthSE(sLength,x(6)))*dS;

            t
            
            % State equations and Lagrangian Formulas
            xdot(1) = x(3);
    
            xdot(2) = x(4);
            
            xdot(3) = -(m2*l1*l2*x(4)^2*sin(x(1)-x(2))*I2COM+m2^2*l1*l2^3*x(4)^2*sin(x(1)-x(2))+m1*g*d1*sin(x(1))*...
                I2COM+m1*g*d1*sin(x(1))*m2*l2^2+m2*g*l1*sin(x(1))*I2COM+m2^2*g*l1*sin(x(1))*l2^2+m2^2*l1^2*l2^2*cos(x(1)-x(2))*...
                x(3)^2*sin(x(1)-x(2))-m2^2*l1*l2^2*cos(x(1)-x(2))*g*sin(x(2)))/(-m2^2*l1^2*l2^2*...
                cos(x(1)-x(2))^2+I1-I2COM+I1 *m2*l2^2+m2*l1^2*I2COM+m2^2*l1^2*l2 ^2) ; %theta1-double-dot
            
            xdot(4) = (l1^2*cos(x(1)-x(2))*m2*l2*x(4)^2*sin(x(1)-x(2))+l1*cos(x(1)-x(2))*m1*g*d1*...
                sin(x(1))+l1^2*cos(x(1)-x(2))*m2*g*sin(x(1))+l1*x(3)^2*sin(x(1)-x(2))*I1+l1^3*...
                x(3)^2*sin(x(1)-x(2))*m2-g*sin(x(2))*I1-g*sin(x(2))*m2*l1^2)*m2*l2/...
                (-m2^2*l1^2*l2^2*cos(x(1)-x(2))^2+I1*I2COM+I1*m2*l2^2+m2*l1^2*I2COM+m2^2*l1^2*l2^2) + (taMoment-sMoment)/I2; %theta 2-double-dot
            
            xdot(5) = TA.getVelocity(aTA,x(5),TA.getNormalizedLengthSE(taLength,x(5)));
            xdot(6) = S.getVelocity(aS,x(6),TA.getNormalizedLengthSE(sLength,x(6)));
            
        end

        function double_pendulum(SP,ivp,duration,fps,TA,S,stimulation)
            % Simulates behaviour of double pendulum with muscle moments
            % included, creates graphs to assist with visualization
            
            % clear All; clf;
            nframes = duration*fps;
            sol = ode45(@(t,x) SP.dynamics(t,x,TA,S,stimulation),[0 duration], ivp);
            t = linspace(0,duration,nframes);
            y = deval(sol,t);
            
            % Define angles, angular velocities, lengths for plotting
            phi1=y(1,:)'; dtphi1=y(3,:)';
            phi2=y(2,:)'; dtphi2=y(4,:)';
            l1=SP.shankLength; l2=SP.ankleFootLength; l3=SP.COMToeLength;
            
            % Coordinates of ankle and COM of foot
            x1 = l1*sin(phi1);
            y1 = -l1*cos(phi1);
            x2 = l1*sin(phi1)+l2*sin(phi2);
            y2 = -l1*cos(phi1)-l2*cos(phi2);
            x3 = l1*sin(phi1)+l2*sin(phi2)+l3*cos(1.43577-phi2);
            y3 = -l1*cos(phi1)-l2*cos(phi2)-l3*sin(1.43577-phi2);
            
            % Plot trajectory of ankle and COM of foot
            figure(1)
            plot(x1,y1,'linewidth',2)
            hold on
            plot(x2,y2,'r','linewidth',2)
            h=gca;
            get(h,'fontSize') 
            set(h,'fontSize',14)
            legend('Ankle','Foot COM')
            xlabel('X (m)','fontSize',14);
            ylabel('Y (m)','fontSize',14);
            title('Ankle and Foot Motion','fontsize',14)
            fh = figure(1);
            set(fh, 'color', 'white'); 
            
            % Plot angle of shank and ankle-foot COM over time
            figure(2)
            plot(t,phi1,'linewidth',2)
            hold on
            plot(t,phi2,'r','linewidth',2)
            h=gca; 
            get(h,'fontSize') 
            set(h,'fontSize',14)
            legend('\theta_1','\theta_2')
            xlabel('Time (s)','fontSize',14);
            ylabel('phi (rad)','fontSize',14);
            title('phi_1 and phi_2','fontsize',14)
            fh = figure(2);
            set(fh, 'color', 'white'); 
            
            % Calculate and plot height relative to first inital height
            figure(3)
            height = zeros(size(x1));
            for i = 1:length(y3)
                height(i) = y3(i) - y3(1);
            end
            plot(t, height)
            xlabel('Time (s)','fontSize',14);
            ylabel('Toe Height (m)','fontSize',14);
            title('Toe Height over Time','fontsize',14)
            axis([0 (SP.duration+0.1) -0.05 0.25]); axis square;

            figure(4)
            h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
            range=1.1*(l1+l2); axis([-range range -2*range range]); axis square;
            set(gca,'nextplot','replacechildren');
            for i=1:length(phi1)-1
                if (ishandle(h)==1)
                    Xcoord=[0,x1(i), x2(i), x3(i)];
                    Ycoord=[0,y1(i), y2(i), y3(i)];
                    set(h,'XData',Xcoord,'YData',Ycoord);
                    drawnow;
                    F(i) = getframe;
                    % pause(t(i+1)-t(i));   
                    pause(0.01);
                end
            end
        end
        
        function isometric(SP,TA,f0M,restingLengthCE,restingLengthSE,stimulation)
            % Simulates isometric force behaviour of the TA
            % included, creates graphs to assist with visualization
            
            L = restingLengthCE + restingLengthSE;
            
            if stimulation
               afun = @(t) SP.getActivationTA(t);
            else               
               afun = @(t) 0*t;
            end
            
            odefun = @(t, y) TA.getVelocity(afun(t), y, TA.getNormalizedLengthSE(L,y));
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45 
            [time, x] = ode45(odefun, [0 SP.duration], 1, OPTIONS);
            
            figure(5)
            subplot(2,1,1)
            plot(time, x*restingLengthCE)
            ylabel('CE length')
            set(gca, 'FontSize', 18)
            xlim([0 (SP.duration+0.1)]);
            subplot(2,1,2)
            plot(time, f0M*TA.forceLengthSE(TA.getNormalizedLengthSE(L, x)));
            xlabel('Time (s)')
            ylabel('Force')
            set(gca, 'FontSize', 18)
            xlim([0 (SP.duration+0.1)]);
        end
    end      
end
