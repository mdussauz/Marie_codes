% Process Lever Signal

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
LeverTrace = TracesOut.Lever{1}; % Lever trace

%% filter the Lever data
% -- adding extra points to allow for the (known) transient of filter at beginning of signal to
% settle
dataToAdd = 15;
StartBuffer = ones(dataToAdd,1)*LeverTrace(1);
LeverTrace = [StartBuffer; LeverTrace];
SampleRate = 500; % Samples/second
filter = 25; %passband frequency of the filter in hertz
FilteredLever = lowpass(LeverTrace,filter,SampleRate);

%% Velocity and acceleration 
LeverVelocity = gradient(LeverTrace);
LeverAcceleration = gradient(LeverVelocity);
VelocityAsIntegralOfAcceleration = cumtrapz(LeverAcceleration); % checking

FilteredLeverVelocity = gradient(LeverFiltered);
FilteredLeverAcceleration = gradient(FilteredLeverVelocity);


%% Plots
start =1;
finish = 20000;

% filtered vs unfiltered
figure(1)
subplot(2,1,1)
plot(LeverTrace(start:finish))
title('Unfiltered lever signal')

subplot(2,1,2)
plot(LeverFiltered(start:finish))
title('Filtered lever signal')

% unfiltered velocity and acceleration
figure(2)
subplot(4,1,1)
plot(LeverTrace(start:finish))
title('Unfiltered lever signal')

subplot(4,1,2)
plot(LeverVelocity(start:finish))
title('Velocity of unfiltered lever signal')

subplot(4,1,3) 
plot(LeverAcceleration(start:finish))
title('Acceleration of unfiltered lever signal')

subplot(4,1,4)   
plot(VelocityAsIntegralOfAcceleration(start:finish))
title('Integral of acceleration = velocity')

%filtered velocity and acceleration
figure(3)
subplot(3,1,1)
plot(FilteredLever(start:finish))
title('Filtered lever signal')

subplot(3,1,2)
plot(FilteredLeverVelocity(start:finish))
title('Velocity of filtered lever signal')

subplot(3,1,3) 
plot(FilteredLeverAcceleration(start:finish))
title('Acceleration of filtered lever signal')

%cutting extra data
figure()
FilteredLeverVelocityCut = FilteredLeverVelocity(dataToAdd+1:end);
plot(FilteredLeverVelocityCut)