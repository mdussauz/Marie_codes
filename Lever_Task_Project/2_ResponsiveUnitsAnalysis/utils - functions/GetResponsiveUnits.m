function [responsive_units, resp_score] = GetResponsiveUnits(SessionName)
% INPUT: 
% SessionName = 'O3/O3_20211005_r0_processed.mat';

% OUTPUT: 
% responsive_units = vector of units where 1 = responsive and 0 = unresponsive
% resp_score = unit x odor x event where 1 = activated and 2 = inhibited

% Function to find cells that are responsive to different events of lever
% task behavior with and without aligning to 1st sniff after trial start: 
% 1) trial on = 0
% 2) odor on = myEvents(:,1)
% 3) trial off = myEvents(:,3)
% 4)reward = most times will be the same as trial off
% 5) 1st TZ entry > 100 ms

% Here neurons are considered responsive if neurons are considered 
% responsive if their average activity in the window after the event 
% exceeded the 99th percentile of the reference response (pre event or iti)
% for all trials.

%% %Settings to be changed depending on analysis performed

% General parameters
BaselineChoice = 'PreEvent'; %PreEvent vs ITI
SortTrials =1; % !!! this is for perturbation and replay trials but once I incorporate those trials it should be deleted from analysis !!!
TZsorted = 1; % to get PSTH average across each TZ type - !!! using all trials doesn't work !!!
ActivityType = 'FR'; %'FR' or 'SpikeCount' - !!! SpikeCount doesn't work !!!

% !!! As of now sniff option is not available for OL and Offset II !!!
sniff_warp = 0; % sniff warp method -% 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency

% Parameters to extract responses:
if strcmp (BaselineChoice, 'PreEvent') % best settings for PreEvent comparison
window = [-350 350]; % time window to extract around event 
threshold = [5 95]; % %[1 99] [3 97] [5 95]; [10 90]; = Thresholds to try 

elseif strcmp (BaselineChoice, 'ITI') %best settings for ITI comparison
window = [-100 100]; % time window to extract around event 
threshold = [5 95]; %[1 99] [3 97] [5 95]; [10 90]; = Thresholds to try 
end 
                             
%% File path

if strcmp(computer,  'PCWIN64')
    ProcessedBehaviorPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
else
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
end
MySession = fullfile(ProcessedBehaviorPath,SessionName);
%% Get the data from processed session loaded
% SingleUnits = spikes aligned to tstart from tstart of this trial TO tstart of next trial trial

[~, ~, handles.TrialInfo, handles.SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, handles.TimestampAdjuster, ~, ~, OpenLoop] = ...
    LoadProcessedDataSession(MySession);  % loads relevant variables

%% Get all spikes, all units aligned to trials - for closed loop and replay (if any)
% AlignedSpikes = % spikes aligned to tstart from previous trial tstop TO next trial t start

% closed-loop:
[handles.AlignedSpikes, handles.Events, handles.whichtetrode, handles.sniffAlignedSpikes, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAlignedSpikeTimes_v2(handles.SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession,'sniffwarpmethod',sniff_warp);

EventsToAlign = 1:5; %1=TrialON ; 2=OdorON ; 3=TrialOFF ; 4=Reward ; 5=1stTZentry>100ms

% replay of any kind:
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay')) 
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(handles.SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    disp('halt flip session')
    EventsToAlign = append(EventsToAlign,6); %6 = perturbation time

    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, handles.SniffAlignedReplaySpikes, handles.SniffAlignedReplayEvents] = ...
        PerturbationReplayAlignedSpikeTimes_v2(handles.SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',sniff_warp);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))
    disp('Offset session')
    EventsToAlign = append(EventsToAlign,6); %6 = perturbation time

    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(handles.SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
end

%% Getting number of units in session 
NbUnit = size(handles.SingleUnits,2);

%% Parameters to analyze responsiveness:
windowsize = length(window(1):window(2)); 
halfwindow = floor(windowsize/2);

baseline_window = [1 halfwindow];
test_window = [halfwindow+1 windowsize];

%% Get SpikeTimes and FR aligned to each events 
%initialize
MeanAlignedActivity = zeros(NbUnit, 3, max(EventsToAlign), windowsize); % unit x odor x event x window length 
AllMeanTrialAlignedFRs = zeros(NbUnit, 3, max(EventsToAlign), windowsize); % unit x odor x event x window length 

for odor = 1:3
    for event = EventsToAlign
        
        for unit = 1:NbUnit
            %  baseline and perturbation trials in closed loop - %dim are trials(or TZ) x time 
            [trialsdone, AlignedFRs, AlignedPerturbationFRs, RawSpikeCounts, RawPerturbationSpikeCounts] = ...
                EventAlignedActivity(unit, odor, handles.AlignedSpikes, handles.Events, handles.TrialInfo, event, window, 'sniffscalar', sniff_warp, 'TZsorted', TZsorted);

            if strcmp(ActivityType, 'SpikeCount')
                AllAlignedActivity.(['odor',num2str(odor)])(:,unit,event,:)= RawSpikeCounts; %dim: trials x unit x event x time
                ThisMeanAlignedActivity = squeeze(mean(RawSpikeCounts, 1)); %take the mean across trials (or TZ)

            elseif strcmp(ActivityType, 'FR')
                AllAlignedActivity.(['odor',num2str(odor)])(:,unit,event,:)= AlignedFRs; %dim: trials x unit x event x time
                ThisMeanAlignedActivity = squeeze(mean(AlignedFRs, 1)); %take the mean across trials (or TZ)
            end
            MeanAlignedActivity (unit,odor,event,:) = ThisMeanAlignedActivity(1:windowsize); %dim: units x odor x event x times
            clearvars ThisMeanAlignedActivity



            %  open loop trials
%             if any(strcmp(handles.TrialInfo.Perturbation,'OL-Replay')) % replay trials
%                 [ActiveAlignedFRs, ActiveRawSpikeCounts, PassiveAlignedFRs, ...
%                     PassiveRawSpikeCounts] = ...
%                     ReplayTrialAlignedActivity(trialsdone, whichunit, whichodor, handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, AlignType, SortTrials);
%             end
%             if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template')) || ...
%                     any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))% passive pertubation trials
% 
%                 [PerturbationReplayAlignedFRs,PerturbationReplayRawSpikeCounts,...
%                     OtherPerturbationReplayAlignedFRs,OtherPerturbationReplayRawSpikeCounts]...
%                     = PerturbationReplayTrialAlignedActivity(trialsdone, whichunit,...
%                     whichodor, handles.AlignedSpikes, Events, TrialInfo, AlignTo, SortTrials, sniffaligned', align_to_sniff, 'sniffscalar',sniff_warp);
% 
%             end
        end
    end
end

%% Window averages of baseline periods and percentile of baseline response probabilities for each unit across trials
% responsive neurons = response probability exceeding the kth (95,97,99) 
% percentile of baseline response probabilities for that neuron;

switch BaselineChoice
    case 'PreEvent' % if we want to compare response to event to activity before event
% AllAlignedFRs.odor(x) = trial x unit x event x time
% baseline_percentile = unit x odor x event
for odor = 1:3
    % get average baseline activity across trials
    % baseline_activity = trial x unit x event
    baseline_activity = squeeze(mean(...
        AllAlignedActivity.(['odor',num2str(odor)])(:,:,:,baseline_window(1,1):baseline_window(1,2)),4,'omitnan')); 
    
    % getting percentile of average baseline resp across trials
    baseline_prctl_low(:,odor,:)   = prctile(baseline_activity,threshold(1,1),1);    % Set the low threshold  
    baseline_prctl_high(:,odor,:)  = prctile(baseline_activity,threshold(1,2),1);    % Set the high threshold 

    % TESTING - get bonus "event" activity - aligned to end of trial ~ animal has settled in TZ
    mean_tz_activity(:,odor) = round(squeeze(mean(baseline_activity(:,:,4), 1)));
end
    case 'ITI' %if we want to compare response to event to activty during ITI (or just before trial on)
for odor = 1:3
    % get average baseline activity just before Trial ON across trials
    % baseline_activity = trial x unit 
    baseline_activity = squeeze(mean(...
        AllAlignedActivity.(['odor',num2str(odor)])(:,:,1,baseline_window(1,1):baseline_window(1,2)),4,'omitnan')); 
    
    % getting percentile of average baseline resp across trials
    baseline_prctl_low(:,odor)   = prctile(baseline_activity,threshold(1,1),1);    % Set the low threshold  
    baseline_prctl_high(:,odor)  = prctile(baseline_activity,threshold(1,2),1);    % Set the high threshold 

    % TESTING - get bonus "event" activity - aligned to end of trial ~ animal has settled in TZ
    mean_tz_activity(:,odor) = squeeze(mean(...
        AllAlignedActivity.(['odor',num2str(odor)])(:,:,4,baseline_window(1,1):baseline_window(1,2)),1,4,'omitnan')); 

end
end

baseline_prctl_low = round(baseline_prctl_low);
baseline_prctl_high = round(baseline_prctl_high);
%% Period data extraction for responsiveness test
% average responses to one odor and one event for that specific window 
% AllMeanCenterAlignedFRs (whichunit,whichodor,AlignTo,:)

mean_event_activity = NaN(NbUnit,3,max(EventsToAlign)); % unit x odor x event 
%calculate mean response for defined time window: 
mean_event_activity(:,:,:) = squeeze(mean(MeanAlignedActivity (:,:,:,test_window(1):test_window(2)),4)); 
mean_event_activity = round(mean_event_activity);

%% Identification of responsive boutons, both enhanced and suppressed for all odors - SINGLE THRESHOLD VALUES 
resp_score = NaN(NbUnit,3, max(EventsToAlign)); % unit x odor x event

% Identification of boutons showing enhanced or suppressed responses 


for odor = 1:3
    for event = EventsToAlign
        for unit = 1:NbUnit

            switch BaselineChoice
                case 'PreEvent'
                    reference_high = baseline_prctl_high(unit,odor,event);
                    reference_low = baseline_prctl_low(unit,odor,event);
                case 'ITI'
                    reference_high = baseline_prctl_high(unit,odor);
                    reference_low = baseline_prctl_low(unit,odor);
            end

            if mean_event_activity(unit,odor,event) > reference_high
                resp_score(unit,odor,event) = 1;
                
            elseif  mean_event_activity(unit,odor, event) < reference_low
                resp_score(unit,odor,event) = 2;
                
            else
                resp_score(unit,odor,event) = 0;
            end
        end
    end
end

% TESTING - Bonus event responsiveness - aligned to end of trial ~ settled in TZ
for odor = 1:3
    for unit = 1: NbUnit
        switch BaselineChoice
            case 'PreEvent'
                reference_high = baseline_prctl_high(unit,odor,5);
                reference_low = baseline_prctl_low(unit,odor,5);
            case 'ITI'
                reference_high = baseline_prctl_high(unit,odor);
                reference_low = baseline_prctl_low(unit,odor);
        end

        if mean_tz_activity(unit,odor) > reference_high
            resp_score(unit,odor,max(EventsToAlign +1)) = 1;

        elseif  mean_tz_activity(unit,odor) < reference_low
            resp_score(unit,odor,max(EventsToAlign +1)) = 2;

        else
            resp_score(unit,odor,max(EventsToAlign +1)) = 0;
        end

    end
end

%% Making response matrix for responsiveness to any event
responsive_units   = NaN(NbUnit,1);

for unit = 1:NbUnit
    
    if any(squeeze(resp_score(unit,:,:))~=0, 'all')
        responsive_units(unit,1) = 1;
    else
        responsive_units(unit,1) = 0;
    end
end