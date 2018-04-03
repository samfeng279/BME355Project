ta_data = csvread('TA_STIM.csv'); 

%TO DO: clip to proper portion of the gait cycle we are modelling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEAST SQUARES POLYNOMIAL REGRESSION CODE TO PRODUCE A SAMPLE DATA SETS
% WITH A GREATER NUMBER OF DATA POINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TO DO: identify the proper gait interval in seconds (duty cycle)
% updating time values to fit with the gait inteval

gait_interval = 2;
for i = 1:size(ta_data,1)
    ta_data(i,1) = ta_data(i,1)*(10^-2)*gait_interval;
end

%plot sample values (curve we are trying to fit)
figure(1)
plot(ta_data(:,1),ta_data(:,2), 'k')
hold on
plot(0,0, 'r') %just to make a red point appear
xlabel('Time (s)')
ylabel('EMG (mV)')
legend('Data', 'Regression')
hold on

% finding the regressed line for each point to map the curve
delta = 5;
interval = 10;
regressions = []; % matrix of time points and their corresponding coeffs
for i = 1:size(ta_data,1)
    % finding start and end points of data interval being regressed
    if i < delta + 1
        s = 1;
        e = i + delta;
    elseif (size(ta_data,1) - i) < delta
        s = i - delta;
        e = size(ta_data,1);
    else
        s = i - delta;
        e = i + delta;
    end
    
    %take care of a zero case, otherwise we only begin when t > 0 (theres
    %probably a cleaner way to do this)
    if i == 1
        x = [0; ta_data(1:9,1)];
        y = [ta_data(1,2); ta_data(1:9,2)];
        coeffs = polyfit(x,y,2);
        regressions = [regressions; 0 coeffs(1) coeffs(2) coeffs(3)];
    end
    
    x = ta_data(s:e,1);
    y = ta_data(s:e,2);
    coeffs = polyfit(x,y,2);
    regressions = [regressions; ta_data(i,1) coeffs(1) coeffs(2) coeffs(3)];
    
    % plot obtained regression
    reg_y = polyval([coeffs(1) coeffs(2) coeffs(3)], x);
    plot(x,reg_y, 'r')
    hold on
end

% using regression to expand sample numbers
% Currently have too few sample points to identify the proper EMG behaviour
% of the muscle, therefore, we need to identify mV at more time points.
% This is done by defining the amount of data points needed in the gait
% interval. Then, it is identified where the point exists (which polynomial
% can be used to identify its activation). These regressed values will tehn
% be passed into the activation algo --> shown step by step below

%creating new time array (can select how many points we want]
new_times = linspace(0, gait_interval, 5000);

simulated_ta = [];
for i = 1:size(new_times,2)
    for j = 1:size(regressions,1)-1
        % finding which interval of the regression data the new time falls into
        if (new_times(i) >= regressions(j,1)) && (new_times(i) <= regressions(j+1,1))
            % now that we've found the interval, we can approx the
            % activation value using either j's or j+1's polynomial
            a = polyval([regressions(j, 2) regressions(j, 3) regressions(j, 4)], new_times(i));
            simulated_ta = [simulated_ta; new_times(i) a];
            break;
        end
    end
end

simulated = simulated_ta(round(0.6*size(simulated_ta,1)):round(1*size(simulated_ta,1)),1:2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACTIVATION MODELLING CODE BASED ON THE SIMULATED SAMPLE DATA SET TO
% PROVIDE MODEL OF HOW TA MUSCLE IS ACTUALLY STIMULATED DURING THE WALKING
% CYCLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulse_width = 0.001; %500-1000 microseconds --> may need to change based on sampling rates
frequency = 50; %20-50Hz

% using frequency to find how many points exist between zeros and peaks
point_dist = (size(new_times,2)/gait_interval)/frequency;
%identifying where peaks and zeros occur, as well as identifying the linear
%relation between adjacents
latest_peak = 1;
latest_zero = 1;
peaks = [];
zeros = [1];
eqns = []; % store info in format: start_index end_index coeff(1) coeff(2)
for i = 1:size(simulated_ta,1)    
    % identify where and when peaks should occur (to create jagged activation)
    if (latest_peak == 1) && ((i - latest_peak) > point_dist/2)
        latest_peak = i;
        peaks = [peaks; i];
    else
        if (i - latest_peak) > point_dist
            latest_peak = i;
            peaks = [peaks; i];
        elseif (i - latest_zero) > point_dist
            latest_zero = i;
            zeros = [zeros; i];
        end  
    end
    
    % we are looking at the ascending curve to peak
    if latest_zero < latest_peak
        % check for empty eqns (check is not necessary for descending case
        % as it is assumed ascending comes first
        if isempty(eqns)
            coeffs = polyfit([simulated_ta(latest_zero,1) simulated_ta(latest_peak,1)], [simulated_ta(latest_zero,2) simulated_ta(latest_peak,2)], 1);
            eqns = [eqns; latest_zero latest_peak coeffs(1) coeffs(2)];
        % checking if the coeffs have already been found for these
        % particular points
        elseif ~any(eqns(:,1) == latest_zero)
            coeffs = polyfit([simulated_ta(latest_zero,1) simulated_ta(latest_peak,1)], [0 simulated_ta(latest_peak,2)], 1);
            eqns = [eqns; latest_zero latest_peak coeffs(1) coeffs(2)];
        end
        
    % we are looking at the descending curve from peak 
    elseif latest_zero > latest_peak
        if ~any(eqns(:,1) == latest_peak)
            coeffs = polyfit([simulated_ta(latest_peak,1) simulated_ta(latest_zero,1)], [simulated_ta(latest_peak,2) 0], 1);
            eqns = [eqns; latest_peak latest_zero coeffs(1) coeffs(2)];
        end
    end
end
     
% finded modelled activation values
activations = [];
for i = 1:size(simulated_ta,1)
    if any(zeros(:) == i)
        % if labelled as a zero, input into the activation matrix with
        % activation of 0
        activations = [activations; simulated_ta(i,1) 0];
    elseif any(peaks(:) == i)
        % if labelled as a peak, input into the activation matrix with
        % sample activation 
        activations = [activations; simulated_ta(i,1) simulated_ta(i,2)];
    else
        % else, based on the location of the point, use the coefficients found 
        % to find the corresponding activation
        % find new value based on linearity between adjacent zero and peak
        for j = 1:size(eqns,1)
            if  i < eqns(j,2) && i > eqns(j,1)
                a = polyval([eqns(j,3) eqns(j,4)],simulated_ta(i,1));
                activations = [activations; simulated_ta(i,1) a];
            end
        end
    end
end


% plot regular sample data and overlay modelled activations
figure(2)
plot(simulated_ta(:,1),simulated_ta(:,2), 'b')
hold on
plot(activations(:,1), activations(:,2), 'r')
xlabel('Time (s)')
ylabel('EMG (smV)')
legend('Simulated', 'Modelled')
coeffs = polyfit(x,y,2);







