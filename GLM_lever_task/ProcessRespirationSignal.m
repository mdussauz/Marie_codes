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
[therm_pks_1,therm_locs_1,therm_w_1,therm_p_1] = findpeaks(ThermistorNormalized, 'MinPeakProminence', 0.01, 'MinPeakDistance', 10);
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
% 
% figure(4)
% subplot(3,1,1)
% plot(ThermistorCut(start:finish))
% 
% subplot(3,1,2)
% plot(ThermistorNormalized(start:finish))
% 
% subplot(3,1,3)
% plot(ThermistorCut(start:finish)); 
% hold on;
% plot(ThermistorNormalized(start:finish))

figure(5)
plot(ThermistorNormalized);
hold on;
plot(therm_locs_1,ThermistorNormalized(therm_locs_1),'o');
hold on;
plot(therm_locs_2,ThermistorNormalized(therm_locs_2),'o');

