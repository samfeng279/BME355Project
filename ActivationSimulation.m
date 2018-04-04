%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATING THE MODEL IN ORDER TO FIND OPTIMUM FREQUENCY FOR MINIMAL
% RELATIVE FATIGUE IN THE TIBIALIS ANTERIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instantiating simulation specific data
data = csvread('TA_STIM.csv'); 
pointNum = 5000;
gaitInterval = 1.41;
swingInterval = [0.6 1];

% finding model simulated data for later use
regressions = getActivationRegression(data, gaitInterval); 
simulated = getSimulatedActivations(regressions, gaitInterval);

% find clipped data for the swing phase
clipped = clipSwingPhase(simulated, swingInterval);
model = getModelledActivations(clipped, 34, gaitInterval, 5000);

figure(1)
plot(clipped(:,1),clipped(:,2), 'b')
hold on
plot(model(:,1), model(:,2), 'r')
xlabel('Time (s)')
ylabel('EMG (smV)')
legend('Simulated', 'Modelled')

% simulate the fatigue for a set of freqencies
fatigue = getFatigue(20, 50, simulated, gaitInterval, pointNum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS FOR SIMULATING FATIGUE AT DIFFERENT FREQUENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
function result = getSimulatedActivations(regressions, gaitInterval)
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
    
    result = simulated;
 end


function result = getModelledActivations(simulatedTA, frequency, gaitInterval, pointNum)
    % using frequency to find how many points exist between zeros and peaks
    delta = round((pointNum/gaitInterval)/frequency);
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

    result = activations;
end 

function result = getSimulated(regressions, gaitInterval, swingInterval)
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
    
    simulated
    
    % clipping data for the swing phase using swingInterval
    min_index = round(swingInterval(1)*size(simulated,1));
    max_index = round(swingInterval(2)*size(simulated,1));
    sim = simulated(min_index:max_index,1:2);
    min_index
    max_index
    sim
    
    result = sim;
end

function result = clipSwingPhase(simulated, swingInterval)
   % clipping data for the swing phase using swingInterval
    minIndex = round(swingInterval(1)*size(simulated,1));
    maxIndex = round(swingInterval(2)*size(simulated,1));
    simulated = simulated(minIndex:maxIndex,1:2);
    result = simulated;
end

function result = getFatigue(minFreq, maxFreq, simulated, gaitInterval, pointNum)
    fatigue = [];
    for i = minFreq:maxFreq
        model = getModelledActivations(simulated, i, gaitInterval, pointNum);
        area = trapz(model(:,1),model(:,2));
        fatigue = [fatigue; i area];
    end
    
    % plot regular sample data and overlay frequency specific modelled activations
    figure(i+1)
    plot(simulated(:,1),simulated(:,2), 'b')
    hold on
    plot(model(:,1), model(:,2), 'r')
    xlabel('Time (s)')
    ylabel('EMG (smV)')
    legend('Simulated', 'Modelled')
    
    % plot frequency vs. fatigue
    figure(i+2)
    plot(fatigue(:,1),fatigue(:,2), 'b')
    xlabel('Frequecy (Hz)')
    ylabel('Relative Fatigue')
    
    result = fatigue;
end






