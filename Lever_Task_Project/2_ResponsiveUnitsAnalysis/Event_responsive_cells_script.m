% Event_responsive_cells_script 
% Script to find cells that are responsive to different events of lever
% task behavior with and without aligning to 1st sniff after trial start: 
% 1) trial on = 0
% 2) odor on = myEvents(:,1)
% 3) trial off = myEvents(:,3)
% 4)reward = most times will be the same as trial off
% 5) 1st TZ entry > 100 ms

%% %Settings to be changed depending on analysis performed

% Parameters to extract responses:
SortTrials =1; % !!! this is for perturbation and replay trials but once I incorporate those trials it should be deleted from analysis !!!

% !!! As of now sniff option is not available for OL and Offset II !!!
sniff_warp = 0; % sniff warp method -% 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency
window = [-500 500]; % time window to extract around event
windowsize = length(window(1):window(2)); 

% Parameters to analyze responsiveness:
baseline = [1 100];
reference = [1 100];
threshold = [1 99]; %[5 95]; [10 90]; = Thresholds to try                                   

%% File path
SessionName = 'O3/O3_20211005_r0_processed.mat';
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

% Getting number of units in session 
NbUnit = size(handles.SingleUnits,2);

%% Get SpikeTimes and FR aligned to each events 

%initialize
MeanAlignedFRs = zeros(NbUnit, 3, 6, windowsize); % unit x odor x event x window length 
AllMeanTrialAlignedFRs = zeros(NbUnit, 3, 6, windowsize); % unit x odor x event x window length 

for whichodor = 1:3
    for AlignTo = EventsToAlign
        
        for whichunit = 1:NbUnit
            %  baseline and perturbation trials in closed loop - %dim are TZ x time 
            [trialsdone, AlignedFRs, AlignedPerturbationFRs, RawSpikeCounts, RawPerturbationSpikeCounts] = ...
                EventAlignedActivity(whichunit, whichodor, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignTo, window, 'sniffscalar', sniff_warp);

            %take the mean across TZ:
            ThisMeanAlignedFR = squeeze(mean(AlignedFRs, 1));
            MeanAlignedFRs (whichunit,whichodor,AlignTo,:) = ThisMeanAlignedFR(1:windowsize); %store in matrix with dimension units x odor x events x times
            clearvars ThisMeanAlignedFR

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

%% Frame window averages of baseline periods. 
% Here I divided the baseline in 4 periods with a window frame size equivalent to the ones I use for the response period 
% MeanAlignedFRs unit x odor x alignto x time
 
baseline_A(:,:) = squeeze(mean(MeanAlignedFRs(:,:,:,baseline(1,1):baseline(1,2)),4,'omitnan')); % 1st window
% baseline_B(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(2,1):reference_baseline(2,2)),3,'omitnan')); % 2nd window
% baseline_C(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(3,1):reference_baseline(3,2)),3,'omitnan')); % 3rd window
% baseline_D(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(4,1):reference_baseline(4,2)),3,'omitnan')); % 4th window

%baseline_all = [baseline_A baseline_B baseline_C baseline_D];              % These values are concatenated in one single matrix  
baseline_all = baseline_A(:,:);

%% Baseline reference values for each unit.
               
baseline_prctl_low_single   = prctile(baseline_all,threshold(1,1),2);    % Set the low threshold  
baseline_prctl_high_single  = prctile(baseline_all,threshold(1,2),2);    % Set the high threshold 

%% Period data extraction for responsiveness test

% % Next I reduce the frame and trial dimensions to one by getting an average value

% what we want is a matrix of dimensions unit x
% tested windows which contains average responses to one odor for that
% specific window 
% AllMeanCenterAlignedFRs (whichodor,whichunit,:)

nCategory = size(reference,1); %how many time windows to test 

resptest_values = NaN(3,6, NbUnit,nCategory); % odor x alignto x unit x category
for nPeriods = 1:size(reference,1)
    %calculate mean response for defined time window: 
    resptest_values(:,:,:,nPeriods) = squeeze(mean(MeanAlignedFRs (:,:,reference(nPeriods,1):reference(nPeriods,2)),3)); 
end 
%% Identification of responsive boutons, both enhanced and suppressed for all odors - SINGLE THRESHOLD VALUES 
resptest_index_single       = NaN(3, 6, NbUnit,nCategory); % odor x align to x unit x time window

% Identification of boutons showing enhanced or suppressed responses for each period
for odor = 1:3
    for alignto = 1:6
    for nCategory = 1:nCategory
        for unit = 1:NbUnit
            if resptest_values(odor, alignto, unit,nCategory) > baseline_prctl_high_single(unit,1)
                resptest_index_single(odor,unit,nCategory) = 1;
                
            elseif  resptest_values(odor, unit,nCategory) < baseline_prctl_low_single(unit,1)
                resptest_index_single(odor,unit,nCategory) = 2;
                
            else
                resptest_index_single(odor, unit,nCategory) = 0;
            end
        end
    end
    end
end


% Classification of boutons by their responsiveness. 
% New array with num Channels x 2 columns. 
% The first column is assignated with a 1 if
% I found at least one enhanced response in one of the periods of each 
% channel row. The second with a 2 if I found at least one suppressed 
% response in each channel row. If no enhanced or suppressed response is 
% found I assign a 0 value in column 1 or 2 respectively. 

resptest_response_single    = NaN(3,NbUnit,2); %odor x unit x enh/inh

for odor = 1:3
    for unit = 1:NbUnit
        
        if find(resptest_index_single(odor,unit,:) == 1,1,'first')
            resptest_response_single(odor,unit,1) = 1;
        else
            resptest_response_single(odor,unit,1) = 0;
        end
        
        if find(resptest_index_single(odor,unit,:) == 2,1,'first')
            resptest_response_single(odor,unit,2) = 2;
        else
            resptest_response_single(odor,unit,2) = 0;
        end
    end
end

%% Making response matrix for responsiveness to any odor
resptest_response_all    = NaN(NbUnit,2); %unit x enh/inh

for unit = 1:NbUnit
    
    if find(resptest_response_single(:,unit,1) == 1,1,'first')
        resptest_response_all(unit,1) = 1;
    else
        resptest_response_all(unit,1) = 0;
    end
    
    if find(resptest_response_single(:,unit,2) == 2,1,'first')
        resptest_response_all(unit,2) = 2;
    else
        resptest_response_all(unit,2) = 0;
    end
end

%% Response class

% new variable named 'resptest_response_class' that is fed with the sum 
% value of each channel row of 'resptest_response_single'. In that way: 
% a enhanced channel corresponds to 1 + 0 = 1; 
% a suppressed channel 0 + 2 = 2; 
% a complex/dual response channel a 1 + 2 = 3; 
% an unresponsive channel a 0 + 0 = 0.

resptest_response_class(:,1) = resptest_response_all(:,1) + resptest_response_all(:,2);
