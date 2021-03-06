%ConstructingOdorIdMove_Test

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

% Lever trace
LeverTrace = TracesOut.Lever{1};

% Motor trace - odor movement
MotorTrace = TracesOut.Motor{1};
MotorTraceAdjusted = TracesOut.Motor{1} + 125; % arbitrary value to remove negative numbers and zeros

% Trials
TrialTrace = TracesOut.Trial{1};
% Trials Sorted by Odor identity 
Odor1 = TrialTrace;
Odor1(Odor1~=1) = 0;
Odor2 = TrialTrace;
Odor2(Odor2~=2) = 0;
Odor2(Odor2==2) = 1;
Odor3 = TrialTrace;
Odor3(Odor3~=3) = 0;
Odor3(Odor3==3) = 1;
Air = TrialTrace;
Air(Air~=4) = 0;

SampleRate = 500; % Samples/second

%% Odor location for each identity
Odor1Location = MotorTraceAdjusted.*Odor1;
Odor2Location = MotorTraceAdjusted.*Odor2;
Odor3Location = MotorTraceAdjusted.*Odor3;
AirLocation = MotorTraceAdjusted.*Air;

%% Plots
start =1;
finish = 20000;
figure();
x = 3;
y = 4;
subplot(x,y,1)
plot(MotorTraceAdjusted(start:finish))
title('Motor Adjusted')

subplot(x,y,5)
plot(Odor1(start:finish))
title('Odor1 Trials')

subplot(x,y,9)
plot(Odor1Location(start:finish))
title('Odor1 Location')

subplot(x,y,6)
plot(Odor2(start:finish))
title('Odor2 Trials')

subplot(x,y,10)
plot(Odor2Location(start:finish))
title('Odor2 Location')

subplot(x,y,7)
plot(Odor3(start:finish))
title('Odor3 Trials')
subplot(x,y,11)
plot(Odor3Location(start:finish))
title('Odor3 Location')        

subplot(x,y,8)
plot(Air(start:finish))
title('Air Trials')
subplot(x,y,12)
plot(AirLocation(start:finish))
title('Air Location')    
