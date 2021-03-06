%ProcessRespirationSignal

%% File path
WhichSession = 'O3_20210918_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3'; %local save
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% Load the relevant variables
load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));

% Get full behavior traces
FirstTrialinAnalysis =2; % ignoring trial 1 as potential misalignment with OpEphys
[TracesOut] = ConcatenateTraces(Traces, FirstTrialinAnalysis:length(TrialInfo.TrialID), SampleRate*startoffset);

SampleRate = 500; % Samples/second
SniffsTrace = TracesOut.Sniffs{1}; % Sniffs

%% filter the thermistor data
% -- adding extra points to allow for the (known) transient of filter at beginning of signal to
% settle
dataToAdd = 15;
StartBuffer = ones(dataToAdd,1)*SniffsTrace(1);
SniffsTrace = [StartBuffer; SniffsTrace];

nqf = SampleRate/2; % Nyquist freq.
[b,a] = butter(3,[0.1 30]/nqf,'bandpass');   % Butterworth filter
ThermistorFiltered = filter(b,a,SniffsTrace);  % filtez

ThermistorFiltered = smoothdata(ThermistorFiltered, 'movmean', 13);

%% Rescale the data
ThermistorCut = ThermistorFiltered(dataToAdd+1:end); % remove extra points
ThermistorNormalized = ThermistorCut - median(ThermistorCut);

%% find points in the thermistor that correspond to the valley
% finding positive peaks
[therm_pks_1,therm_locs_1,therm_w_1,therm_p_1] = findpeaks(ThermistorNormalized, 'MinPeakProminence', 0.01, 'MinPeakDistance', 10);
% finding negative peaks
[therm_pks_2,therm_locs_2,therm_w_2,therm_p_2] = findpeaks(-ThermistorNormalized, 'MinPeakProminence', 0.01, 'MinPeakDistance', 10);

%% Plots
start =1;
finish = 25000;

% figure(1) 
% subplot(2,1,1)
% plot(SniffsTrace(start:finish))
% title('Unfiltered thermistor signal')
% 
% subplot(2,1,2)
% plot(ThermistorFiltered(start:finish))
% title('Filtered thermistor signal')
% 
% figure(2)
% plot(ThermistorCut)
% title('Filtered thermistor signal after cutting data buffer')
% 
% figure(3)
% plot(ThermistorNormalized)
% title('Filtered and normalized thermistor signal')
% 
% figure(4)
% subplot(3,1,1)
% plot(ThermistorCut(start:finish))
% title('Thermistor signal without extra datapoints')
% 
% subplot(3,1,2)
% plot(ThermistorNormalized(start:finish))
% title('Thermistor signal without extra datapoints')
% 
% subplot(3,1,3)
% plot(ThermistorCut(start:finish)); 
% hold on;
% plot(ThermistorNormalized(start:finish))
% title('Both plotted on top of each other')

figure(5)
plot(ThermistorNormalized);
hold on;
plot(therm_locs_1,ThermistorNormalized(therm_locs_1),'o');
hold on;
plot(therm_locs_2,ThermistorNormalized(therm_locs_2),'o');
title('Positive and negative peaks of thermistor signal')

%% combine both positive and negative peaks 
therm_pks_combined = sort([therm_pks_1(:)',therm_pks_2(:)']);
therm_locs_combined = sort([therm_locs_1(:)',therm_locs_2(:)']);

% -- Plot
figure(6)
plot(ThermistorNormalized);
hold on;
plot(therm_locs_combined,ThermistorNormalized(therm_locs_combined),'o');
title('Both negative and positive peaks')

%% Interpolating breathing signal based on peaks of inhalation and exhalation
% -- Adding 1st and final value of Thermistor signal otherwise
% interpolation is wrong at the beginning and at the end
x = [1, therm_locs_combined, length(ThermistorNormalized)];
y = [ThermistorNormalized(1); ThermistorNormalized(therm_locs_combined); ThermistorNormalized(end)];

% -- Defining the stretch of data
xx = 1:1:length(ThermistorNormalized); 

% -- Trying different kind of interpolating functions
yy = spline(x,y,xx);
p = pchip(x,y,xx);
m = makima(x,y,xx);

% Plots
figure(7)
plot(x,y,'o',xx,yy)
title('spline function results')
figure(8)
plot(x,y,'o',xx,p)
title('pchip function results')
figure(9)
plot(x,y,'o',xx,m)
title('makima function results')
%pchip function seems to work best for what I am trying to achieve

%% keeping signal between 2 values to remove slow oscillations

%-- defining positive peaks as equal to 2 and negative peaks as equal to 1
y_locs_1 = ones(length(therm_locs_1),1)*2;
y_locs_2 = ones(length(therm_locs_2),1);

%-- getting the index of sorted combined peaks
[B,I] = sort([therm_locs_1(:)',therm_locs_2(:)']);

%-- combining the new 1/2 peaks and sorting them based on index I
y_combined = [y_locs_1(:)', y_locs_2(:)'];
y_combined = y_combined(I);

%-- adding start and end of signal for better interpolation after
%scaling between 1 and 2
scaled_y = normalize(y, 'range', [1 2]);
y_combined = [scaled_y(1); y_combined(:); scaled_y(end)];

%-- interpolating between peaks to get an oscillating signal between peaks
interpolated_thermistor = pchip(x,y_combined,xx);

% Plots
figure(10)
plot(x,y_combined, 'o',xx, interpolated_thermistor);
hold on;
plot(therm_locs_combined,ThermistorNormalized(therm_locs_combined),'o', xx,ThermistorNormalized);
legend('New peaks','Interpolated signal','Old peaks','Original signal')

