%BinningMotorTrace

%% what session

WhichSession = 'O3_20211005_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% Load relevant variables

load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.SampleRate = SampleRate;

%% Load Template and Open Loop trial info 

x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
u = unique(TrialInfo.Perturbation(x));
handles.PerturbationList.String = u{1};
for y = 2:size(u,1)
    handles.PerturbationList.String = [handles.PerturbationList.String,'; ',u{y}];
end

% get start and stop TS of the template
templateStart = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'first'));
templateEnd   = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'last'));

% get start and stop of the replays
replays = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Replay')));

%% Get full behavior traces for Baseline Closed Loop 
FirstTrialinAnalysis =2; % we chose to ignore trial 1 as potential misalignment with OpEphys
LastClosedLoopBaseline = replays(1)-1;
[TracesOut] = ConcatenateTraces(Traces, FirstTrialinAnalysis:LastClosedLoopBaseline, SampleRate*startoffset);

%% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

%% Create New variable combining odor and motor move
% Odor Trials 
TrialTrace = TracesOut.Trial{1};
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
Air(Air==4) = 1;
% --
MotorTrace = TracesOut.Motor{1}; % Motor trace - odor movement

% -- Binning 
Binsize = 50;
BehaviorBinsize = Binsize/(1000/SampleRate);
myMotor = reshape(MotorTrace,BehaviorBinsize,[]);

% -- Adjust Motor to remove zeros
MotorTraceAdjusted = MotorTrace + 125; % arbitrary value to remove negative numbers and zeros

subplot(2,1,1)
plot(MotorTrace(1:20000))
subplot(2,1,2)
plot(myMotor(1:20000))

