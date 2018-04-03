% for loop which goes through various frequencies
	% perform peak finding ~
	% finding area under stimulus
	% for loop which iterates through all peaks
		% val = multiply period by amplitude of stimulus (bxh/2*2)
		% integrate += val
	
	% store frequency and final integration in 2D array

data = csvread('TA_STIM.csv'); 
pointNum = 5000;
gaitInterval = 2;
swingInterval = 1;
minFreq = 20;
maxFreq = 50;
fatigue = [];
for i = minFreq:maxFreq
    %getting the modelled data for the specified frequency
    regressions = getActivationRegression(data, gaitInterval);
    simulatedTA = getSimulatedActivations(regressions, gaitInterval, swingInterval);
    model = getModelledActivations(simulatedTA, i, gaitInterval, pointNum);
    
    delta = (pointNum/gaitInterval)/i; %finding the curret delta to identify the proper peaks
    % re-finding to data peaks --> this section could be optimized
    latestPeak = 1; latestZero = 1;
    areas = []; % holder array for each half-isoceles triangle area
    for j = 1:size(model,1)    
        % identify where and when peaks should occur (to create jagged activation)
        if (latestPeak == 1) && ((j - latestPeak) > delta/2)
            latestPeak = j;
        else
            if (j - latestPeak) > delta
                latestPeak = j;
            elseif (j - latestZero) > delta
                latestZero = j;
            end  
        end

        % Indentifying the zero-peak relationship to find the area under
        % each linear relationship
        if latestZero < latestPeak
            % check for empty eqns (check is not necessary for descending case
            % as it is assumed ascending comes first
            if isempty(areas)
                height = model(latestPeak,2);
                base = abs(model(latestPeak,1) - model(latestZero,1));
                area = base*height/2;
                areas = [areas; area];
            % checking if the coeffs have already been found for these
            % particular points
            elseif ~any(areas(:,1) == latestZero)
                height = model(latestPeak,2);
                base = abs(model(latestPeak,1) - model(latestZero,1));
                area = base*height/2;
                areas = [areas; area];
            end
        % looking at the descending curve from peak 
        elseif latestZero > latestPeak
            if ~any(areas(:,1) == latestPeak)
                height = model(latestPeak,2);
                base = abs(model(latestZero,1) - model(latestPeak,1));
                area = base*height/2;
                areas = [areas; area];
            end
        end
    end
    
    fatigue = [fatigue; i sum(areas)];
    
    % plot regular sample data and overlay modelled activations
    figure(i)
    plot(simulatedTA(:,1),simulatedTA(:,2), 'b')
    hold on
    plot(model(:,1), model(:,2), 'r')
    xlabel('Time (s)')
    ylabel('EMG (smV)')
    legend('Simulated', 'Modelled')
end

% plot freq vs. fatigue
figure(i+1)
plot(fatigue(:,1),fatigue(:,2), 'b')
xlabel('Frequecy (Hz)')
ylabel('Relative Fatigue')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS
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

    % TO DO: clip data for the swing phase using swingInterval
    % (will also need to alter code in below function)

    result = simulated;
end

function result = getModelledActivations(simulatedTA, frequency, gaitInterval, pointNum)
    % using frequency to find how many points exist between zeros and peaks
    delta = (pointNum/gaitInterval)/frequency;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATING THE MODEL IN ORDER TO FIND OPTIMUM FREQUENCY FOR MINIMAL
% RELATIVE FATIGUE IN THE TIBIALIS ANTERIOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





