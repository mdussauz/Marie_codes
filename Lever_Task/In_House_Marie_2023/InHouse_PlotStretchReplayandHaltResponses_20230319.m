% Get processed file path
WhichSession = 'O8\O8_20220704_r0_processed.mat';
MyUnits = [];
[Paths] = WhichComputer();
handles.WhereSession.String = fullfile(Paths.Local.Behavior_processed,WhichSession);


%% get the data loaded
MySession = handles.WhereSession.String;
[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables



% %% Get all spikes, all units aligned to trials
% [handles.AlignedSpikes, handles.Events, handles.whichtetrode] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
%     size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);
% 
% if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay'))
%     [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
%         ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
%         ReplayTTLs,handles.TrialInfo,handles.Events);
% end
% 
% if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
%     [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
%         PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
%         ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
% end
% 
% handles.NumUnits.String = num2str(size(SingleUnits,2));
% handles.CurrentUnit.Data(1) = 1;

%% get the processed data loaded
% load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
%                 'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
%                 'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
% 
% OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% which Units to use
N = size(SingleUnits,2); % #units
if isempty(MyUnits)
    MyUnits = 1:N;
end

%% Get the replay traces and spikes
[MyTraces,timestamps,PSTH,Raster] = ...
ProcessHaltReplayTrialsMD(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
                        'plotfigures', 1, 'plotephys', 1, ...
                        'PlotOpenLoop', 1, ...
                            'whichunits', MyUnits);