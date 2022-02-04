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

%% Plots
start =1;
finish = 25000;

figure()
subplot(2,1,1)
plot(SniffsTrace(start:finish))
title('Unfiltered thermistor signal')

subplot(2,1,2)
plot(ThermistorFiltered(start:finish))
title('Filtered thermistor signal')

%%
figure()
ThermistorCut = ThermistorFiltered(dataToAdd+1:end);
plot(ThermistorCut)