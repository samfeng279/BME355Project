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
            
            odefun = @(t, x) HillTypeMuscle.getVelocity(afun(t), x, m.getNormalizedLengthSE(L, x));
            
            OPTIONS = odeset('AbsTol', 1e-6, 'RelTol', 1e-5);  % use this as the final argument to ode45 
            [time, x] = ode45(odefun, [0 2], 1, OPTIONS);
            
            figure
            subplot(2,1,1)
            plot(time, x*m.restingLengthCE)
            ylabel('CE length')
            set(gca, 'FontSize', 18)
            subplot(2,1,2)
            
            % TO DO: should go up at the end!
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
            
            alpha = 0; % pennation angle
            beta = 0.1; % damping coefficient (see damped model in Millard et al.)
            
            seFL = HillTypeMuscle.forceLengthSE(lT);
            ceFL = HillTypeMuscle.forceLengthCE(lM);
            peFL = HillTypeMuscle.forceLengthPE(lM);
            
            fun = @(vM) (a*ceFL*HillTypeMuscle.forceVelocityCE(vM) + peFL + beta*vM)*cos(alpha) - seFL;
            result = fzero(fun, 0); 
        end
        
        function result = forceLengthCE(lM)
            % Normalized force-length curve of contractile element. 
            % 
            % lM: contracile element length
            % result: force-length scale factor
            
            result = HillTypeMuscle.forceLengthRegression.eval(lM);            
        end
        
        function result = forceVelocityCE(vM)
            % Normalized force-velocity curve of contractile element.  
            % 
            % vM: contracile element velocity
            % result: force-velocity scale factor
            
            result = max(0, HillTypeMuscle.forceVelocityRegression.eval(vM));
        end
        
        function result = forceLengthSE(lT)
            % Normalized force-length curve of series elastic element. 
            % 
            % lT: normalized length of tendon (series elastic element)
            % result: force produced by tendon
            
            sT = 1;
            result = zeros(size(lT));
            for i = 1:size(lT)
                if (lT(i) >= sT) 
                    result(i) = 10 * (lT(i) - sT) + 240 * (lT(i) - sT) ^ 2;
                else
                    result(i) = 0;
                end
            end
        end
        
        function result = forceLengthPE(lM)
            % Normalized force-length curve of parallel elastic element.
            % 
            % lM: normalized length of contractile element
            % result: force produced by parallel elastic element
                      
            sM = 1;
            
            result = zeros(size(lM));
            for i = 1:size(lM,2)
                if (lM(i) >= sM) 
                    result(i) = 3 * (lM(i) - sM) ^ 2 / (0.6 + lM(i) - sM);
                else
                    result(i) = 0;
                end
            end
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
    %   from Winters et al. (2011) Fig. 3C, normalized so that max force is
    %   ~1 and length at max force is ~1. 
    
    
    % TO DO: fix up so that the values don't go negative!
    
    data = [
        40.34655403660477, 36.87001724439568
        45.71414290355635, 68.30080223927234
        56.87268304467716, 100.74225876590111
        57.79604961720776, 92.32499437682753
        53.347995234882006, 89.85579686601858
        67.38709919276235, 52.602904056181615
        72.84599171935788, 26.485392247519542
        75.28253317671758, 13.079749081548357
        62.63264439057306, 80.55581935870845
        43.76827531052409, 57.32736860520333
        40.31856313364823, 17.72423962212281
        37.324036354851344, 9.467922925049322
        50.777663925890764, 91.7491815160073
        54.65040528494907, 100.70427111188866
        65.71347645348597, 67.84495039112284
        42.09431934621248, 32.341489016069715
        39.188180508001565, 24.54252367980405
        38.296636926332276, 14.726713818009102
        38.26739642949375, 14.726213980456293
        39.30422612650889, 3.917726738809847
        41.611059738918186, 1.7919176267712373
        41.24992710702355, 14.777197410841438
        41.25092678212914, 15.460975183065443
        40.29532068744325, 21.826406417914143
        41.39188097201743, 31.87364106665325
        41.3845500212431, 26.859270737010434
        42.842242937711916, 23.921225601679396
        43.368571880805405, 23.93022267762973
        43.717125267621356, 22.34073925974053
        43.383900232424466, 34.41481518506487
        43.384899907530055, 35.09859295728887
        42.51801497846534, 42.14930147701995
        43.39922858404352, 44.89940769249989
        42.81941702280093, 48.30829980256408
        42.81608477244896, 46.02904056181734
        43.20254250701855, 50.36613100742255
        44.39540482676464, 46.283957713743064
        45.50754338173428, 46.98672931297335
        45.506210481593484, 46.07502561667451
        45.44439723756446, 43.79476670082221
        46.6738310049234, 44.72746357433829
        45.45939236414833, 54.05143328418262
        45.16765384583344, 54.50228675680398
        43.47187164171645, 54.587259140779224
        44.29976924166313, 60.86921750431111
        46.37867692999776, 62.842076325194284
        47.700663950882635, 67.08119861045162
        47.343779938186756, 62.972533926473886
        49.011321320570815, 63.570839477169955
        46.719149609710186, 75.72538924849425
        46.42474529111372, 74.35283532851821
        46.42191287831456, 72.41546497388356
        45.95256541623972, 71.38180091470275
        47.444497205075024, 71.86314447804466
        47.16575446313282, 81.20310898957837
        47.48839960346221, 81.89238497488316
        48.134356334191395, 83.72678879364207
        48.897608277309885, 85.79111788668678
        48.95142412049418, 82.60115462474698
        49.009238664100835, 82.14630245170315
        49.711177200743094, 82.27226151500759
        49.46942243770775, 76.91200359883038
        50.63637651096727, 75.10858970834482
        50.875632086238646, 78.75940319396199
        53.63223618990496, 84.27661010171693
        53.56675747048877, 79.48916602104327
        51.05474054265697, 81.26958738410013
        51.11022251101726, 79.21925374253368
        50.65137163755114, 85.36525629170518
        50.70985263122819, 85.36625596681074
        50.80073975957815, 87.53305175817854
        50.304817600946365, 88.3222952540425
        49.65986054532278, 87.17166920750759
        49.7751564075009, 86.03403893734537
        51.71302659968844, 91.53725039362206
        51.594398487158344, 90.39562142303748
        51.594398487158344, 90.39562142303748
        52.64672314831013, 90.18568965086348
        55.75587933921476, 96.84852422962541
        56.04828430760004, 96.85352260515339
        53.93730370962771, 92.94279359208258
        53.645565191312826, 93.393647064704
        54.34900324061347, 94.54527278634447
        53.472788010563235, 95.2140554319846
        53.563508526395594, 97.2668882613151
        53.886153666724994, 97.95616424661984
        54.232707703329766, 94.99912528428263
        56.37526137338699, 100.50583560342892
        57.51564075008956, 100.52532926798796
        57.74923149976259, 100.30140204433559
        57.27838452502938, 98.2420713268188
        56.987978906855275, 99.6046284957389
        58.50481926707154, 97.12343488366287
        58.37985987887271, 91.65121335565942
        59.3155557777056, 91.66720815734885
        59.384200134956146, 98.619948516732
        58.91851814826849, 100.09346962237288
        60.02432542756938, 96.46564866418413
        59.61495847182999, 96.45865093844503
        59.8505485717142, 97.60227925924079
        60.61513341497347, 100.57831204858422
        60.576895842184626, 94.42381226101514
        61.368555219553656, 95.91882638142607
        61.42303751280835, 93.18471496763561
        61.416539624622004, 88.74015944817938
        61.29424603670475, 85.0913453127734
        61.40454352335492, 80.53482618149104
        61.43011854480628, 78.02814085422239
        62.42271261839903, 96.96248719166277
        62.17879189263489, 90.12071076900008
        63.69963095327353, 90.37462824582002
        63.13939636284876, 87.17416839527152
        63.63581835903332, 86.72681378551971
        63.81292746524046, 87.86944243120985
        64.51469938936513, 87.88143853247692
        64.5073684385908, 82.86706820283405
        63.802930714184555, 81.03166470896957
        63.5098592957289, 80.57081448529232
        63.394896658585964, 81.93637067952909
        63.91456110097553, 77.38684927398589
        66.36976316030625, 76.74505785619675
        65.14199551812328, 76.9519906030539
        66.07202659135781, 73.09324469547397
        65.60484509201177, 73.54109914277853
        65.71547580369715, 69.21250593557085
        65.3876656753222, 64.99037812710873
        63.304259449012406, 59.94051933121739
        63.29476253550929, 53.444630495089086
        64.81676788377112, 54.49628870617056
        64.81460192104234, 53.014770199685074
        65.68565216304702, 48.81313573088744
        66.9640700105799, 43.25094344338089
        67.24614500287407, 36.19023817259392
        67.75397995651414, 63.5493464623998
        68.3919392863986, 59.91352810336639
        67.11135547613694, 63.994201884387564
        66.82294920817402, 66.72431459775567
        70.48225993218871, 49.692849823807194
        69.35979140112798, 41.92437457826196
        69.36079107623357, 42.60815235048602
        71.39729587883939, 35.57743733286674
        70.45393580419697, 30.319146277459765
        70.22001182948875, 30.315147577037408
        68.46158331875475, 27.550046234973536
        72.2013678887695, 25.562692125059357
        73.26802122643474, 35.15357508809632
        74.99812560917702, 18.544972883812704
        74.3471705029199, 13.291680203933595
        73.35299361040997, 13.274685727138547
        73.24536192404136, 19.654612251018307
        73.24336257383018, 18.287056706570297
        76.32919301227102, 8.99507660010488
    ];
    
    [optForce, index] = max(data(:,2));
    optLength = data(index,1);
    normalizer = [optLength, optForce];
    data = bsxfun(@rdivide, data, normalizer);
    
    length = data(:,1);
    force = data(:,2);

    centres = linspace(0.6, 1.2, 6);
    width = .3;
    lambda = .01;
    result = Regression(length, force, centres, width, lambda, 1);

end
