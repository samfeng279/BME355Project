classdef HillTypeMuscle < handle
    % Damped Hill-type muscle model adapted from Millard et al. (2013). The
    % dynamic model is defined in terms of normalized length and velocity. 
    % To model a particular muscle, scale factors are needed for force, CE
    % length, and SE length. These are given as constructor arguments. 
    
    properties (Access = public)
        f0M;               
        restingLengthCE;    
        restingLengthSE;    
    end
    
    properties (Constant)
        % curve fits for CE force-length and force-velocity 
        forceLengthRegression = getCEForceLengthRegression();
        forceVelocityRegression = getCEForceVelocityRegression();
    end
    
    methods (Access = public)
        function m = HillTypeMuscle(f0M, restingLengthCE, restingLengthSE)
            % The arguments are scale factors model is normalized 
            % f0M: maximum isometric force
            % restingLengthCE: actual length of CE (m) that corresponds to 
            %   normalized length of 1
            % restingLengthSE: % actual length of SE (m) that corresponds 
            %   to normalized length of 1
            
            m.f0M = f0M;
            m.restingLengthCE = restingLengthCE;
            m.restingLengthSE = restingLengthSE;
        end
        
        function result = getNormalizedLengthSE(m, muscleTendonLength, normalizedLengthCE)
            % Calculates the normalized length of the series elastic
            % element. 
            % 
            % muscleTendonLength: physical (non-normalized) length of the
            %   full muscle-tendon complex (typically found from joint 
            %   angles and musculoskeletal geometry)
            % normalizedLengthCE: normalized length of the contractile
            %   element (the state variable of the muscle model)
            % result: normalized length of the series elastic element
            
            result = (muscleTendonLength - m.restingLengthCE*normalizedLengthCE) / m.restingLengthSE;
        end
        
        function result = getForce(m, length, normalizedLengthCE)
            % length: muscle-tendon length (m)
            % normalizedLengthCE: normalized length of CE (the state
            %   variable)
            result = m.f0M * m.forceLengthSE(m.getNormalizedLengthSE(length, normalizedLengthCE));
        end
                
        function simulateIsometric(m)
            L = m.restingLengthCE + m.restingLengthSE;
            afun = @(t) t>0.5;
            
            odefun = @(t, y) HillTypeMuscle.getVelocity(afun(t), y, getNormalizedLengthSE(m,L,y));
            
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45 
            [time, x] = ode45(odefun, [0 2], 1, OPTIONS);
            
            figure
            subplot(2,1,1)
            plot(time, x*m.restingLengthCE)
            ylabel('CE length')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(time, m.f0M*HillTypeMuscle.forceLengthSE(m.getNormalizedLengthSE(L, x)));
            xlabel('Time (s)')
            ylabel('Force')
            set(gca, 'FontSize', 18)
        end
        
    end
    
    methods (Static)
        
        function result = getVelocity(a, lM, lT)
            % Calculates normalized velocity of contractile element. 
            % 
            % a: activation (between 0 and 1)
            % lM: normalized length of contractile element (this is
            %   \tilde{l}^M in the paper and the lecture slides)
            % lT: normalized length of series elastic element (this is
            %   \tilde{l}^T in the paper and the lecture slides)
            % result: normalized velocity
            
            beta = 0.1; % damping coefficient (see damped model in Millard et al.)
            fL = HillTypeMuscle.forceLengthCE(lM);
            fPE = HillTypeMuscle.forceLengthPE(lM);
            fT = HillTypeMuscle.forceLengthSE(lT);
           
            alpha = 0; % given pennation angle
            % Ignore f0M as it is redundant within function
            input = @(vM) (a*fL*HillTypeMuscle.forceVelocityCE(vM)+fPE+beta*vM)*cos(alpha)-fT;
            result = fzero(input, 0);

        end
        
        function result = forceLengthCE(lM)
            % Normalized force-length curve of contractile element. 
            % lM: contracile element length
            % result: force-length scale factor
            
            result = HillTypeMuscle.forceLengthRegression.eval(lM);            
        end
        
        function result = forceVelocityCE(vM)
            % Normalized force-velocity curve of contractile element.  
            % vM: contracile element velocity
            % result: force-velocity scale factor
            
            result = max(0, HillTypeMuscle.forceVelocityRegression.eval(vM));
        end
        
        function result = forceLengthSE(lT)
            % Normalized force-length curve of series elastic element. 
            % lT: normalized length of tendon (series elastic element)
            % result: force produced by tendon
            indices = lT <1;
            result = arrayfun(@(lT) 10*(lT-1) + 240*((lT - 1).^2), lT);
            result(indices) = 0;
        end
        
        function result = forceLengthPE(lM)
            % Normalized force-length curve of parallel elastic element.
            % lM: normalized length of contractile element
            % result: force produced by parallel elastic element
            indices = lM <1;
            result = arrayfun(@(lM) 3*((lM-1).^2)/(0.6+lM-1), lM);
            result(indices) = 0;
        end
        
        function plotCurves()
            % Plot force-length, force-velocity, SE, and PE curves. 
            
            lM = 0:.01:1.8;
            vM = -1.2:.01:1.2;
            lT = 0:.01:1.07;
            figure
            subplot(2,1,1), hold on
            plot(lM, HillTypeMuscle.forceLengthCE(lM), 'r')
            plot(lM, HillTypeMuscle.forceLengthPE(lM), 'g')
            plot(lT, HillTypeMuscle.forceLengthSE(lT), 'b')
            legend('CE', 'PE', 'SE', 'location', 'northwest')
            xlabel('Normalized length')
            ylabel('Force scale factor')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            plot(vM, HillTypeMuscle.forceVelocityCE(vM), 'k')
            set(gca, 'FontSize', 18)
            xlabel('Normalized CE velocity')
            ylabel('Force scale factor')
        end
        
    end
    
end

function result = getCEForceVelocityRegression()
    % result: regression model of contractile element force-length curve
    %   based on data in Millard et al.
    
    data = [-1.0028395556708567, 0.0024834319945283845
    -0.8858611825192801, 0.03218792009622429
    -0.5176245843258415, 0.15771090304473967
    -0.5232565269687035, 0.16930496922242444
    -0.29749770052593094, 0.2899790099290114
    -0.2828848376217543, 0.3545364496120378
    -0.1801231103040022, 0.3892195938775034
    -0.08494610976156225, 0.5927831890757294
    -0.10185137142991896, 0.6259097662790973
    -0.0326643239546236, 0.7682365981934388
    -0.020787245583830716, 0.8526638522676352
    0.0028442725407418212, 0.9999952831301149
    0.014617579774061973, 1.0662107025777694
    0.04058866536166583, 1.124136223202283
    0.026390887007381902, 1.132426122025424
    0.021070257776939272, 1.1986556920827338
    0.05844673474682183, 1.2582274002971627
    0.09900238201929201, 1.3757434966156459
    0.1020023112662436, 1.4022310794556732
    0.10055894908138963, 1.1489210160137733
    0.1946227683309354, 1.1571212943090965
    0.3313459588217258, 1.152041225442796
    0.5510200231126625, 1.204839508502158];

    velocity = data(:,1);
    force = data(:,2);

    centres = linspace(-.5, 0, 6);
    width = .3;
    lambda = .01;
    result = Regression(velocity, force, centres, width, lambda, 1);
end

function result = getCEForceLengthRegression()
    % result: Regression model of normalized CE force-length curve. Data 
    %   from Winters et al. (2011) Fig. 2C, normalized so that max force is
    %   ~1 and length at max force is ~1. 
    
    data = [41.83374083, 2
            39.48655257, 3.714285714
            37.50611247, 9.428571429
            38.38630807, 14.85714286
            41.46699267, 14.57142857
            41.39364303, 15.71428571
            40.4400978, 17.14285714
            40.4400978,	21.42857143
            39.3398533,	24.28571429
            43.88753056, 22
            43.44743276, 23.42857143
            42.93398533, 23.42857143
            41.46699267, 26.85714286
            41.39364303, 31.42857143
            42.12713936, 32
            43.44743276, 34.85714286
            40.4400978, 36.85714286
            42.49388753, 41.71428571
            43.44743276, 44.57142857
            44.47432763, 45.42857143
            45.57457213, 43.42857143
            45.57457213, 45.42857143
            45.57457213, 46.57142857
            46.74816626, 44.57142857
            42.8606357, 46.28571429
            42.8606357, 48.28571429
            43.22738386, 50.57142857
            43.44743276, 54.28571429
            45.42787286, 53.42857143
            45.28117359, 54
            43.81418093, 56.85714286
            44.47432763, 60.57142857
            46.45476773, 62.57142857
            47.40831296, 62.57142857
            49.02200489, 62.85714286
            47.77506112, 66.57142857
            45.79462103, 67.71428571
            45.94132029, 70.57142857
            46.45476773, 71.42857143
            47.40831296, 71.42857143
            46.45476773, 73.71428571
            46.67481663, 75.42857143
            50.63569682, 74.28571429
            49.46210269, 76.28571429
            50.78239609, 78.28571429
            51.07579462, 78.57142857
            51.00244499, 80.57142857
            49.68215159, 82
            49.02200489, 81.71428571
            47.55501222, 81.71428571
            47.04156479, 80.57142857
            48.21515892, 83.42857143
            48.94865526, 85.42857143
            49.82885086, 85.14285714
            49.60880196, 86.57142857
            50.70904645, 85.14285714
            50.70904645, 87.14285714
            50.26894866, 87.71428571
            50.78239609, 90.57142857
            51.58924205, 89.71428571
            51.80929095, 91.14285714
            52.61613692, 89.42857143
            53.27628362, 88.85714286
            53.56968215, 83.42857143
            53.56968215, 78.85714286
            53.86308068, 92.28571429
            53.56968215, 92.28571429
            53.49633252, 94.85714286
            53.56968215, 96.57142857
            53.86308068, 96.85714286
            54.22982885, 94.28571429
            54.59657702, 99.25827815
            55.77017115, 96.28571429
            56.13691932, 96.28571429
            56.43031785, 99.68211921
            56.79706601, 99.89403974
            57.31051345, 99.04635762
            57.53056235, 98.83443709
            57.82396088, 99.42857143
            58.85085575, 99.14285714
            58.48410758, 96.28571429
            57.75061125, 91.14285714
            58.33740831, 90.85714286
            59.36430318, 91.14285714
            59.58435208, 95.71428571
            59.36430318, 98
            59.87775061, 96.28571429
            60.53789731, 99.42857143
            60.53789731, 93.71428571
            61.41809291, 94.85714286
            62.29828851, 96.28571429
            61.34474328, 92
            62.22493888, 89.42857143
            61.41809291, 88.28571429
            61.27139364, 84.28571429
            63.10513447, 86.28571429
            63.69193154, 85.71428571
            63.69193154, 87.14285714
            63.61858191, 90
            64.42542787, 87.14285714
            64.49877751, 82
            63.83863081, 80.28571429
            63.47188264, 80
            62.59168704, 80
            61.41809291, 80
            61.41809291, 77.42857143
            63.91198044, 76.28571429
            65.08557457, 76
            66.40586797, 75.71428571
            66.11246944, 72
            65.52567237, 72.57142857
            65.74572127, 66.57142857
            65.67237164, 68.28571429
            65.30562347, 64.28571429
            66.77261614, 66.28571429
            67.1393643, 63.14285714
            67.72616137, 62.57142857
            68.38630807, 60
            63.32518337, 59.42857143
            63.32518337, 52.85714286
            64.79217604, 52
            64.79217604, 54
            67.43276284, 51.71428571
            65.67237164, 48.28571429
            66.99266504, 43.14285714
            70.51344743, 48.57142857
            69.3398533, 41.71428571
            71.3202934, 34.57142857
            67.3594132, 35.71428571
            68.4596577, 27.14285714
            70.07334963, 29.42857143
            70.29339853, 29.42857143
            73.3007335, 34.57142857
            72.93398533, 25.71428571
            72.34718826, 24.57142857
            73.3007335, 18.57142857
            73.3007335, 17.42857143
            75.13447433, 17.71428571
            73.37408313, 12.57142857
            74.32762836, 12.57142857
            75.35452323, 12.28571429
            76.38141809, 8.571428571
            57.09046455, 99.71428571
            59.80440098, 96.85714286
            63.32518337, 80.85714286];
        
    % normalize the data set
    [Mf, I] = max(data(:,2)); %maximum force and its index from data
    Ml = data(I,1); %length at that index
    n_data(:,2) = data(:,2) ./ Mf; %divide all forces by maximum force
    n_data(:,1) = data(:,1) ./ Ml; %divide all lengths by length @index
    
    length = n_data(:,1); %bookkeeping
    force = n_data(:,2); %bookkeeping
    
    %Center of normalized data
    %[Mfn, In] = max(n_data(:,2)); %maximum force and its index from normalized data
    %Mln = n_data(In,1); %normalized length at that index
    %centres = [Mln, Mfn]; %centres of normalized length and force
    centres = [1,1]; %force and length normalized about 1, 1
    
    %Call Regression.m file
    width = .3;
    lambda = .01;
    result = Regression(length, force, centres, width, lambda, 1);

end
