%MeasureDistanceFromFakeNostrils
% For info distance between mice nostrils is roughly 1.8mm
% Our Target Zone is 3 mm
% The odor can be moved at -/+50  mm from center 
% Approximately 120 motoer loc on either side so 0.4mm between 2 points
% So 10 points is 9 dist between 2 points times 0.4 = 3.6mm

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
% MotorTraceAdjusted = TracesOut.Motor{1} + 125; % arbitrary value to remove negative numbers and zeros
% % --
% Odor1Location = MotorTraceAdjusted.*Odor1;
% Odor2Location = MotorTraceAdjusted.*Odor2;
% Odor3Location = MotorTraceAdjusted.*Odor3;
% AirLocation = MotorTraceAdjusted.*Air;

%% create lign for each faske nostril and calculate distance
x_length = length(MotorTrace);
rightNostril = ones(x_length,1)*(-5); % need to check whether negative values are on the right/left side
leftNostril = ones(x_length,1)*5;

odorDist_rightNostril = abs(MotorTrace - rightNostril);
odorDist_leftNostril = abs(MotorTrace - leftNostril);

odor1Dist_right = odorDist_rightNostril .*Odor1;
odor2Dist_right = odorDist_rightNostril .*Odor2;
odor3Dist_right = odorDist_rightNostril .*Odor3;

odor1Dist_left = odorDist_leftNostril .*Odor1;
odor2Dist_left = odorDist_leftNostril .*Odor2;
odor3Dist_left = odorDist_leftNostril .*Odor3;

%% Plots
figure(1)
plot(rightNostril(1:20000)); hold on; plot(leftNostril(1:20000)); hold on; 
plot(MotorTrace(1:20000))
figure(2)
plot(odorDist_rightNostril(1:20000)); hold on; plot(odorDist_leftNostril(1:20000)); hold on; 

figure(3)
plot(odor1Dist_right(1:20000)); hold on; plot(odor2Dist_right(1:20000)); hold on; plot(odor3Dist_right(1:20000)); 

%% Creating tine shifted array 
odorDistRightShifted{1} = odorDist_rightNostril;
odorDistLeftShifted{2} = odorDist_leftNostril;

for i = 1:25 
    odorDistRightShifted{i+1} = circshift(odorDist_rightNostril,i,1);
    odorDistLeftShifted{i+1} = circshift(odorDist_leftNostril,i,1);
end 
clear odorDist_rightNostril; clear odorDist_leftNostril;
odorDist_rightNostril = cat(2,odorDistRightShifted{:});
odorDist_leftNostril = cat(2,odorDistLeftShifted{:});

%% Plots
figure(4)
plot(odorDist_rightNostril(1:20000,1)); hold on ;plot(odorDist_rightNostril(1:20000,2));

%%
odor1Dist_right = odorDist_rightNostril .*Odor1;
odor2Dist_right = odorDist_rightNostril .*Odor2;
odor3Dist_right = odorDist_rightNostril .*Odor3;

odor1Dist_left = odorDist_leftNostril .*Odor1;
odor2Dist_left = odorDist_leftNostril .*Odor2;
odor3Dist_left = odorDist_leftNostril .*Odor3;

figure(3)
plot(odor1Dist_right(1:20000, 1)); hold on; plot(odor1Dist_right(1:20000,2));