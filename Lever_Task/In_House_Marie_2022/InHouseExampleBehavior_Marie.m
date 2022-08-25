SessionPath = 'O5/O5_20210927_r1_processed.mat';
datapath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
            
[~, ~, TrialInfo] = LoadProcessedDataSession(MySession); % to get target zone entry time

%% Trial Aligned spikeTimes
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

%% Behavior plotting
for i = 1:12
    whichTZ = i;
    
    % control trials
    controlTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(~strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')));
    controlTrials(:,3) = Events(controlTrials,5);

    for j = 1:numel(controlTrials(:,1))
        Traces.HaltFlip(controlTrials(j,1)) = {0*Traces.Trial{controlTrials(j,1)}};
        Traces.TargetZone(controlTrials(j,1)) = {TrialInfo.TargetZoneType(controlTrials(j,1)) + ...
            0*Traces.Trial{controlTrials(j,1)}};
    end
        
    % perturbation trials
    perturbationTrials = intersect(find(TrialInfo.TargetZoneType==whichTZ),...
                              find(strcmp(TrialInfo.Perturbation(:,1),'Halt-Flip')));
    perturbationTrials(:,3) = Events(perturbationTrials,5);
    
    for j = 1:numel(perturbationTrials(:,1))
        foo = 0*Traces.Trial{perturbationTrials(j,1)};
        idx = TrialInfo.Perturbation{perturbationTrials(j,1),2};
        idx(:,1:2) = idx(:,1:2) + SampleRate*startoffset;
        foo(idx(1):idx(2)) = -idx(3);
        Traces.HaltFlip(perturbationTrials(j,1)) = {foo};
        Traces.TargetZone(perturbationTrials(j,1)) = ...
                {TrialInfo.TargetZoneType(perturbationTrials(j,1)) + ...
                0*foo};
    end
end

B1 = figure;
figure(B1);
[TracesOut] = ConcatenateTraces(Traces, 1:size(Traces.Lever,2), SampleRate*startoffset);

timestamps = (1:size(TracesOut.Lever{1},1))'/SampleRate;
TZ = [TracesOut.TargetZone{1} TracesOut.TargetZone{1}];
TZ = [TargetZones(TZ(:,1),1) TargetZones(TZ(:,1),3)];
Trial = TracesOut.Trial{1};
%Trial(Trial~=whichOdor) = 0;
PlotBehaviorPerturbations(timestamps,TracesOut.Lever{1},...
    TracesOut.Sniffs{1},TracesOut.Licks{1},TracesOut.Rewards{1},...
    Trial,...
    TZ, ...
    TracesOut.HaltFlip{1},5);

% make odor ON time available on the plot
traceoffset = ((TrialInfo.TraceIndices(1,1))-1)/SampleRate;
plot(-traceoffset+(TrialInfo.OdorStart(:,1)+TrialInfo.SessionTimestamps(:,1)),-0.5,'.k');      
set(gcf,'Position',[0 100 1000 300]);