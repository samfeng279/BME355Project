classdef SwingPhaseModel < handle
    % Simple model of standing postural stability, consisting of foot and
    % body segments, and two muscles that create moments about the ankles,
    % tibialis anterior and soleus. 
    
    properties (Access = public)
        % TO DO: initial values
        
        shankAngle = 7*pi/4; % Initial angle of shank from pi/2
        shankAngVel = pi/12;  % Initial angular velocity of shank
        footAngle = 2*pi+0.6; % Initial angle of COM of foot from pi/2
        footAngVel = 0; % Initial angular velocity of COM of foot
        shankMass = 3.285; 
        footMass = 1.095; 
        shankLength = 0.4318; 
        ankleFootLength = 0.07136; % Distance from ankle to COM of foot
        duration = 0.6; % Duration of simulation
        fps = 100; % Frames per sec - changes resolution of animation  
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
        
        function aS = getActivationS(theta)
            % TO DO: function for soleus activation energy
            aS = 0;
        end
        
        function aTA = getActivationTA(theta)
            % TO DO: function for tibialis activation energy
            aTA = 900;
        end
    end
    
    methods (Access = public) 
        function simulate(SP)
            clc; figure;
            
            % TO DO: establish max iso forces, initial CE lengths of each
            % muscle, resting lengths of muscles
            
            % Initialize TA and soleus muscles
            restLengthTA = SP.tibialisLength(pi/2);
            TA = HillTypeMuscle(222000, .6*restLengthTA, .4*restLengthTA);
            restLengthS = SP.soleusLength(pi/2);
            S = HillTypeMuscle(8000, .6*restLengthS, .4*restLengthS);
            
            % Set initial parameters for ODE
            ivp=[SP.shankAngle; SP.footAngle; SP.shankAngVel; SP.footAngVel; 1; 1];
            
            % Simulate double pendulum behaviour
            SP.double_pendulum(ivp, SP.duration, SP.fps, TA, S);
        end 
        
        function xdot = dynamics(SP,t,x,TA,S)
            % x: [shankAngle; footAngle; shankAngVel; footAngVel; CE TA length; CE soleus length];
            
            xdot = zeros(6,1);
            m1 = SP.shankMass; 
            m2 = SP.footMass; 
            l1 = SP.shankLength; 
            d1 = 0.494 * l1;
            l2 = SP.ankleFootLength;
            g = 9.81; 
%             I1 = m1 * d1 ^ 2;
%             I2 = m2 * l2 ^ 2;
            I1 = 100;
            I2 = 80;
            
            % TO DO: decide what to pass into these activation functions (can be angle/time?)
            aTA = SwingPhaseModel.getActivationTA(t);
            aS = SwingPhaseModel.getActivationS(t);
            
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

            
            
            % TO DO:::::::::fix equations so mass of shank is in the MIDDLE
            % of the shank rather than at the end
            xdot(1) = x(3);
        
%             xdot(3)=-((g*(2*m1+m2)*sin(x(1))+m2*(g*sin(x(1)-2*x(2))+2*(l2*x(4)^2+...
%                 l1*x(3)^2*cos(x(1)-x(2)))*sin(x(1)-x(2))))/...
%                 (2*l1*(m1+m2-m2*cos(x(1)-x(2))^2)));
%             xdot(4)=(((m1+m2)*(l1*x(3)^2+g*cos(x(1)))+l2*m2*x(4)^2*cos(x(1)-x(2)))*...
%                 sin(x(1)-x(2))+taMoment-sMoment)/(l2*(m1+m2-m2*cos(x(1)-x(3))^2));
    
            xdot(2) = x(4);
            
            xdot(3) = -(m2*l1*l2*x(4)^2*sin(x(1)-x(2))*I2+m2^2*l1*l2^3*x(4)^2*sin(x(1)-x(2))+m1*g*d1*sin(x(1))*...
                I2+m1*g*d1*sin(x(1))*m2*l2^2+m2*g*l1*sin(x(1))*I2+m2^2*g*l1*sin(x(1))*l2^2+m2^2*l1^2*l2^2*cos(x(1)-x(2))*...
                x(3)^2*sin(x(1)-x(2))-m2^2*l1*l2^2*cos(x(1)-x(2))*g*sin(x(2)))/(-m2^2*l1^2*l2^2*...
                cos(x(1)-x(2))^2+I1-I2+I1 *m2*l2^2+m2*l1^2*I2+m2^2*l1^2*l2 ^2) ; %theta1-double-dot
            
            xdot(4) = (l1^2*cos(x(1)-x(2))*m2*l2*x(4)^2*sin(x(1)-x(2))+l1*cos(x(1)-x(2))*m1*g*d1*...
                sin(x(1))+l1^2*cos(x(1)-x(2))*m2*g*sin(x(1))+l1*x(3)^2*sin(x(1)-x(2))*I1+l1^3*...
                x(3)^2*sin(x(1)-x(2))*m2-g*sin(x(2))*I1-g*sin(x(2))*m2*l1^2)*m2*l2/...
                (-m2^2*l1^2*l2^2*cos(x(1)-x(2))^2+I1*I2+I1*m2*l2^2+m2*l1^2*I2+m2^2*l1^2*l2^2) + (taMoment-sMoment)/I2 ; %theta 2-double-dot
            
            xdot(5) = TA.getVelocity(aTA,x(5),TA.getNormalizedLengthSE(taLength,x(5)));
            xdot(6) = S.getVelocity(aS,x(6),TA.getNormalizedLengthSE(sLength,x(6)));
            
        end

        function double_pendulum(SP,ivp,duration,fps,TA,S)
            % Simulates behaviour of double pendulum with muscle moments
            % included, creates graphs to assist with visualization
            
            clear All; clf;

            nframes = duration*fps;
            sol = ode45(@(t,x) SP.dynamics(t,x,TA,S),[0 duration], ivp);
            t = linspace(0,duration,nframes);
            y = deval(sol,t);
                        
            phi1=y(1,:)'; dtphi1=y(3,:)';
            phi2=y(2,:)'; dtphi2=y(4,:)';
            l1=SP.shankLength; l2=SP.ankleFootLength;
            
            % Coordinates of ankle and COM of foot
            x1 = l1*sin(phi1);
            y1 = -l1*cos(phi1);
            x2 = l1*sin(phi1)+l2*sin(phi2);
            y2 = -l1*cos(phi1)-l2*cos(phi2);
            
            % Plot trajectory of ankle and COM of foot
            figure(1)
            plot(x1,y1,'linewidth',2)
            hold on
            plot(x2,y2,'r','linewidth',2)
            h=gca; 
            get(h,'fontSize') 
            set(h,'fontSize',14)
            xlabel('X','fontSize',14);
            ylabel('Y','fontSize',14);
            title('Chaotic Double Pendulum','fontsize',14)
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
            xlabel('time','fontSize',14);
            ylabel('phi','fontSize',14);
            title('phi_1 and phi_2','fontsize',14)
            fh = figure(2);
            set(fh, 'color', 'white'); 
            
            % Plot angle of shank and ankle-foot COM over time
            figure(3)
            height = zeros(size(x1));
            height(1) = x1(1) + cos(phi2(1))*0.2386;
            for i = 1:length(x1)
               if cos(phi2) <= 1
                   height(i) = (x1(i) + cos(phi2(i))*0.2386);
               else
                   height(i) = (x1(i) - cos(phi2(i))*0.2386);
               end
            end
            plot(t, height)

            figure(4)
            h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
            range=1.1*(l1+l2); axis([-range range -range range]); axis square;
            set(gca,'nextplot','replacechildren');
            for i=1:length(phi1)-1
                if (ishandle(h)==1)
                    Xcoord=[0,x1(i), x2(i)];
                    Ycoord=[0,-l1*cos(phi1(i)),-l1*cos(phi1(i))-l2*cos(phi2(i))];
                    set(h,'XData',Xcoord,'YData',Ycoord);
                    drawnow;
                    F(i) = getframe;
                    % pause(t(i+1)-t(i));   
                    pause(0.01);
                end
            end
        end
    end      
    
end

