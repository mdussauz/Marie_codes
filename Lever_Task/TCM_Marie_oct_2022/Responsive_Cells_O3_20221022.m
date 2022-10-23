%Responsive_Cells_O3_20221022

%% Select a processed file to analyze

% [WhichSession, SessionPath] = uigetfile(); % to be used later for more flexibility

WhichSession = 'O3_20211005_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);


%% get the data loaded
MySession = handles.WhereSession.String;
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession);  % loads relevant variables

% SingleUnits contains the spiketimes of all units from start of trial up
% to the start of next trial 

% startoffset = Time window (in seconds) preceding trial start used for extracting Traces
startoffset = 1;

%% Get all spikes, all units aligned to trials - both for closed-loop and replays

% closed-loop
[handles.AlignedSpikes, handles.Events, handles.whichtetrode] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

% trialalignedspikes = spiketimes relative to the most recent trial''s'' start timestamp');

% Perturbations - Might be useful at some point but not used right now
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
end

handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.CurrentUnit.Data(1) = 1;

%% Getting number of units in session > will be useful to normalize results .i.e responsiveness 
handles.NumUnits.String = num2str(size(SingleUnits,2));

N = size(SingleUnits,2);
MyUnits = (1:N);

%% Trial order
%Information I need here 
whichUnit = 6; %  SPECIFY 

whichodor = 3; % SPECIFY

for whichodor = 1:3
    for whichUnit = 1:N
        AlignTo = 5; % to get aligned to TZ entry 
        [AlignedFRs, AlignedRawSpikeCounts] = AlignedFRToEvent(whichUnit, whichodor, AlignTo);
        %take the mean across TZ
        %store it in matrix that has dimension odor x units x times 
        
        AlignTo = 1; %Trial ON so I can get FR during ITI
        % repeat store in baseline matrix
    end
end

%%
%%%%
%% FUNCTIONS 
function [meanAlignedFRs, meanAlignedRawSpikeCounts] = AlignedFRToEvent(whichUnit, whichodor, AlignTo)

thisUnitSpikes = handles.AlignedSpikes(:,whichUnit);
% get the trial sorting order
whichTrials = intersect(find(cellfun(@isempty, handles.TrialInfo.Perturbation(:,1))), ...
    find(handles.TrialInfo.Odor==whichodor));
whichTrials = [whichTrials handles.TrialInfo.TargetZoneType(whichTrials) handles.TrialInfo.Duration(whichTrials)];
whichTrials = sortrows(whichTrials,2);

for tz = 1:12
    whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
end

% also collect perturbation trials
perturbationTrials = intersect(find(~cellfun(@isempty, handles.TrialInfo.Perturbation)), ...
    find(handles.TrialInfo.Odor==whichodor));
perturbationTrials = intersect(find(~strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay')), ...
    perturbationTrials);
perturbationTrials = [perturbationTrials handles.TrialInfo.TargetZoneType(perturbationTrials) handles.TrialInfo.Duration(perturbationTrials)];
perturbationTrials = sortrows(perturbationTrials,2);
for tz = 1:12
    perturbationTrials(perturbationTrials(:,2)==tz,:) = sortrows(perturbationTrials(perturbationTrials(:,2)==tz,:),3);
end

allTrials = vertcat(whichTrials, perturbationTrials);
%% USER - select with event to align to

myEvents = handles.Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,1);
    case 2 % odor ON
        odorON = myEvents(:,1);
        myEvents(:,1) = 0; % replace odorON with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - odorON;
        Xlims = [-1.2 -1];
        Offset = odorON;
    case 3 % trial OFF
        TrialOFF = myEvents(:,3);
        myEvents(:,3) = 0; % replace TrialOFF with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - TrialOFF;
        Xlims = [-1.2 -1] - 4;
        Offset = TrialOFF;
    case 4 % reward
        Reward = myEvents(:,3);
        myEvents(:,2) = 0; % replace Reward with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - Reward;
        Xlims = [-1.2 -1] - 4;
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        Offset = myEvents(:,4);
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
        %         case 6 % perturbation start
        %             Offset = myEvents(:,5);
        %             % offset all events with ON timestamp
        %             myEvents = myEvents - Offset;
        %             Xlims = [-1.2 -1];
end

%% Calculate PSTHs for specific event
AlignedFRs = []; RawSpikeCounts = [];
BinOffset = Xlims(1)*1000; % would be window start here

for TZ = 1:12
    thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1)); %spiketimes for all trials of 1 tz type
    Events2Align = Offset(find(whichTrials(:,2)==TZ),1); %time at which event is happening for all trials of tz type
    [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500); %resulting psth is an average of all trials of tz type
    AlignedFRs(TZ,1:numel(myFR)) = myFR;
    RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
end

%entries_done = TZ; % might be needed when looking at perturbations

% calculate PSTH for perturbed trials - commented for now
% AlignedPerturbationFRs = [];
% for TZ = 1:12
%     thisTZspikes = thisUnitSpikes(perturbationTrials(find(perturbationTrials(:,2)==TZ),1));
%     Events2Align = Offset(x+find(perturbationTrials(:,2)==TZ),1);
%     [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
%     AlignedPerturbationFRs(TZ,1:numel(myFR)) = myFR;
%     RawSpikeCounts(TZ+entries_done,1:numel(myPSTH)) = myPSTH;
% end

%x = size(allTrials,1);

end
    
