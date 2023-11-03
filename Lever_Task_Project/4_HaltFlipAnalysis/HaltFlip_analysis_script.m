%HaltFlip_analysis_script

%% Sessions
%SessionName = 'O3/O3_20210929_r0_processed.mat';
SessionName = 'O8/O8_20220704_r0_processed.mat';
%SessionName = 'O9/O9_20220702_r0_processed.mat';
%SessionName = 'S1/S1_20230327_r0_processed.mat';
%SessionName = 'S3/S3_20230327_r0_processed.mat'; %%not sorted
%SessionName = 'S6/S6_20230710_r0_processed.mat';
%SessionName ='S7/S7_20230608_r0_processed.mat';
%SessionName ='S11/S11_20230801_r0_processed.mat';
%SessionName ='S12/S12_20230731_r0_processed.mat';

%% Path
if strcmp(computer,  'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionName);

%% get the processed data loaded

[TracesOut, ColNames, handles.TrialInfo, handles.SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop, handles.Tuning] = ...
    LoadProcessedDataSession(MySession); 

Nb_unit = size(handles.SingleUnits,2);
ChosenUnits = 1:Nb_unit;

%% get the closed loop tuning curve

[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionName, ChosenUnits, 'tuningbins', 15);

%% check that its a halt session - if not disable Halt only options
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')) || ...
            any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        handles.OnlyHaltRelated.Value = 1; 
        handles.OnlyHaltRelated.Enable = 'on';
        handles.OdorList = mode(...
            [ handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')); ...
            handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))] ...
            );
else
    handles.OnlyHaltRelated.Value = 0; 
    handles.OnlyHaltRelated.Enable = 'off'; 
    handles.OdorList = [1 2 3];
end

%% Closed Loop: Get all spikes, all units aligned to trials  

[handles.AlignedSniffs, handles.sniffAlignedSpikes, handles.trialAlignedSpikes, ...
    handles.whichtetrode, handles.Events, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAndSniffAlignedSpikeTimes(handles.SingleUnits,TTLs,size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

%% same for replays
if ~isempty(OpenLoop)
    
    % sniffs
    [handles.ReplayAlignedSniffs, handles.SniffAlignedReplaySpikes, handles.ReplayInfo] = ...
        SniffAlignedSpikeTimes_Replays(handles.SingleUnits,TTLs,ReplayTTLs,handles.TrialInfo,OpenLoop,MySession);
    
    % regular trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayTrialInfo] = ...
            PerturbationReplayAlignedSpikeTimes_v2(handles.SingleUnits,TTLs,...
            ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',0);
    end
else
    handles.ReplayAlignedSniffs = [];
end

%% also for passive tuning
handles.TuningSniffs = PassiveTuningSniffs(handles.Tuning,MySession);

%% pseudorandomtuning trials
if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
    [handles.PseudoRandomTuningSpikes] = ...
        TrialAlignedSpikeTimes_Tuning(handles.SingleUnits,handles.Tuning.TTLs);
end

%%
for whichUnit = 1:Nb_unit
for i = 1:numel(handles.OdorList)
    whichodor = handles.OdorList(i);
    % add the trial aligned plot
    AlignType = 6; % perturbation start
    myXlim = [-1.2 6];
    
    % plot baseline trials
    [trialsdone, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = PlotFullSession(whichUnit, whichodor, handles.trialAlignedSpikes, handles.Events, ...
        handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotspikes', 0, ...
        'trialfilter', handles.PlotSelectTrials.Value);
    
    % passive halts
    if ~isempty(handles.ReplayAlignedSniffs)
        [perturbationreplaysadded, PassiveReplayFRs, PerturbationReplayFRs, BinOffset] = AddPerturbationReplay2FullSession_v2(trialsdone, whichUnit, whichodor, handles.ReplayAlignedSpikes, ...
            handles.ReplayEvents, handles.ReplayTrialInfo, handles.ReplayTrialInfo.InZone, AlignType, handles.SortReplay.Value, ...
            'trialfilter', handles.PlotSelectTrials.Value, 'plotspikes', 0);
        
        trialsdone = trialsdone + perturbationreplaysadded;
    end
    
    % add tuning trials
    if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
        AlignType = 1000 + haltlocation;
        LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
        [trialsdone, FR, BinOffset] = PlotRandomTuningTrials(trialsdone, whichUnit, whichodor, handles.PseudoRandomTuningSpikes, ...
            handles.TuningTiming, handles.Tuning.extras.sequence, AlignType, LocationDuration, myXlim, 'plotspikes', 0);
    else
        [trialsdone, FR, BinOffset] = PlotTuningTrials(trialsdone, whichUnit, whichodor, handles.SingleUnits, handles.Tuning.TTLs, ...
            'plotspikes', 0, 'selectlocation', haltlocation);
    end

end
end