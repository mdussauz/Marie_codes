% Event_responsive_cells_script 
% Script to find cells that are responsive to different events of lever
% task behavior with and without aligning to 1st sniff after trial start: 
% 1) trial on = 0
% 2) odor on = myEvents(:,1)
% 3) trial off = myEvents(:,3)
% 4)reward = most times will be the same as trial off
% 5) 1st TZ entry > 100 ms

%% File path
SessionName = 'O3/O3_20211005_r0_processed.mat';
if strcmp(computer,  'PCWIN64')
    ProcessedBehaviorPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
else
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
end
MySession = fullfile(ProcessedBehaviorPath,SessionName);

%Settings to be changed depending on analysis performed
sniff_warp = 0; % sniff warp method -% 0 - don't align, 1 - simple warp, 2 - complex warp, 3 - fixed latency
if sniff_warp~=0
    align_to_sniff = 1;
else
    align_to_sniff = 0;
end 


%% get the data from processed session loaded
%MySession = handles.WhereSession.String;
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession);  % loads relevant variables

handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjuster;

%% Get all spikes, all units aligned to trials - for closed loop and replay (if any)

% closed-loop:
[handles.AlignedSpikes, handles.Events, handles.whichtetrode, handles.sniffAlignedSpikes, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAlignedSpikeTimes_v2(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession,'sniffwarpmethod',sniff_warp);

% replay of any kind:
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo, handles.SniffAlignedReplaySpikes, handles.SniffAlignedReplayEvents] = ...
        PerturbationReplayAlignedSpikeTimes_v2(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',sniff_warp);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Offset-II-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
end

% Getting number of units in session - useful to normalize results .i.e % of responsiveness: 
N = size(SingleUnits,2);

%% Get SpikeTimes aligned to each events 

SortTrials =1;
window_length = 200; 

AllMeanCenterAlignedFRs = zeros(3, 6, N, window_length); % odor x event x unit x window length 
AllMeanTrialAlignedFRs = zeros(3, 6, N, window_length); % odor x event xunit x window length 

for whichodor = 1:3
    for AlignTo = 1:5 %%% !!! Need to change back to 6 and correct that perturbation inclusion problem !!!
        for whichunit = 1:N

            %  baseline and perturbation trials in closed loop - %dim are TZ x time 
            [trialsdone, AlignedFRs, AlignedPerturbationFRs, RawSpikeCounts, RawPerturbationSpikeCounts] = ...
                TrialAlignedActivity(whichunit, whichodor, handles.AlignedSpikes, handles.Events, handles.TrialInfo, AlignTo, 'sniffaligned', align_to_sniff, 'sniffscalar', sniff_warp);

            if length(AlignedFRs) < window_length % in case the vector is too short to be stored, add zeros to make it the right length
                AlignedFRs(:,window_length) = 0;
            end
            %store each tz responses individually:
            AllCenterAlignedFRs(whichunit,whichodor,AlignTo, 1:12,1:window_length) = AlignedFRs(:,1:window_length);
            %take the mean across TZ:
            MeanCenterAlignedFRs = squeeze(mean(AlignedFRs, 1));
            AllMeanCenterAlignedFRs (whichodor,AlignTo,whichunit,:) = MeanCenterAlignedFRs(1:window_length); %store it in matrix that has dimension odor x units x times
            clearvars MeanCenterAlignedFRs

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

%% Defining parameters to analyze responsiveness
% input dataset for Diego has dim time x unit x trial
num_Channel                 = N; 
%bin                         = 25;                                          % Set binning for baseline distribution plots 
% target zone has been entered 100 ms before (at bin 1100 here) 
%reference                   = [600 700; 700 800; 800 900; 900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
%trying reduced windows - gives the same result:
%reference                   = [800 900; 900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
% trying even smaller window - final settings - 600ms (300 bins) around target zone entry on each side: 
%reference                   = [900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
reference                   = [1 100];
% bin 600 is when the trial is ON 
%reference_baseline          = [100 200; 200 300; 300 400; 400 500];        % You have to divide your baseline in frame windows equivalent to the ones you want to evaluate.
reference_baseline          = [1 100];
reference_threshold         = [1 99; 5 95; 10 90];                         % Here you can feed possible thresholds you want to try. To get an idea of which threshold to pick, use 'Responsive_bouton_threshold_test_20221021.m'                                   
threshold                   = 1;                                           % Pick the desired threshold by defining the used row of reference threshold.

%% Frame window averages of baseline periods. 
% Here I divided the baseline in 4 periods with a window frame size equivalent to the ones I use for the response period 
AllMeanCenterAlignedFRs = reshape(AllMeanCenterAlignedFRs, [N,3,6,window_length]); %unit x odor x alignto x time
 
baseline_A(:,:) = squeeze(mean(AllMeanCenterAlignedFRs(:,:,:,reference_baseline(1,1):reference_baseline(1,2)),4,'omitnan')); % 1st window
% baseline_B(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(2,1):reference_baseline(2,2)),3,'omitnan')); % 2nd window
% baseline_C(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(3,1):reference_baseline(3,2)),3,'omitnan')); % 3rd window
% baseline_D(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(4,1):reference_baseline(4,2)),3,'omitnan')); % 4th window

%baseline_all = [baseline_A baseline_B baseline_C baseline_D];              % These values are concatenated in one single matrix  
baseline_all = baseline_A(:,:);
%baseline_all = baseline_all'; % to correct for how I organize things vs Diego - dim are channel x (tz x odor x 4 window)
%% Baseline reference values for each channel.
               
set_threshold                   = reference_threshold(threshold,:); 

baseline_Roi_prctl_low_single   = prctile(baseline_all,set_threshold(1,1),2);    % Set the low threshold  
baseline_Roi_prctl_high_single  = prctile(baseline_all,set_threshold(1,2),2);    % Set the high threshold 

baseline_Roi_mean               = mean(baseline_all,2,'omitnan');

baseline_Roi_median             = median(baseline_all,2,'omitnan');

%% Period data extraction for responsiveness test
% Currently I am dividing the after cue period in 5 parts to test if in one
% of those I found any changes. You must pick proper frame ranges by
% modifing the 'reference' array and include or remove more dataset windows
% depending on your data.

%window_size = reference(1,2)-reference(1,1);% I shouldn't need the +1 - % here for me would be 100 bins = 200ms 

% % Next I reduce the frame and trial dimensions to one by getting an average value

% what we want is a matrix of dimensions unit x
% tested windows which contains average responses to one odor for that
% specific window 
% AllMeanCenterAlignedFRs (whichodor,whichunit,:)

nCategory = size(reference,1); %how many time windows to test 

resptest_values = NaN(3,6, num_Channel,nCategory); % odor x alignto x unit x category
for nPeriods = 1:size(reference,1)
    %calculate mean response for defined time window: 
    resptest_values(:,:,:,nPeriods) = squeeze(mean(AllMeanCenterAlignedFRs (:,:,reference(nPeriods,1):reference(nPeriods,2)),3)); 
end 
%% Identification of responsive boutons, both enhanced and suppressed for all odors - SINGLE THRESHOLD VALUES


num_Category = nCategory; % the number of tested windows  
resptest_index_single       = NaN(3, 6, num_Channel,num_Category); % odor x align to x unit x time window

% Identification of boutons showing enhanced or suppressed responses for each period
for odor = 1:3
    for alignto = 1:6
    for nCategory = 1:num_Category
        for nChannel = 1:num_Channel
            if resptest_values(odor, alignto, nChannel,nCategory) > baseline_Roi_prctl_high_single(nChannel,1)
                resptest_index_single(odor,nChannel,nCategory) = 1;
                
            elseif  resptest_values(odor, nChannel,nCategory) < baseline_Roi_prctl_low_single(nChannel,1)
                resptest_index_single(odor,nChannel,nCategory) = 2;
                
            else
                resptest_index_single(odor, nChannel,nCategory) = 0;
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

resptest_response_single    = NaN(3,num_Channel,2); %odor x unit x enh/inh

for odor = 1:3
    for nChannel = 1:num_Channel
        
        if find(resptest_index_single(odor,nChannel,:) == 1,1,'first')
            resptest_response_single(odor,nChannel,1) = 1;
        else
            resptest_response_single(odor,nChannel,1) = 0;
        end
        
        if find(resptest_index_single(odor,nChannel,:) == 2,1,'first')
            resptest_response_single(odor,nChannel,2) = 2;
        else
            resptest_response_single(odor,nChannel,2) = 0;
        end
    end
end

%% Making response matrix for responsiveness to any odor
resptest_response_all    = NaN(num_Channel,2); %unit x enh/inh

for nChannel = 1:num_Channel
    
    if find(resptest_response_single(:,nChannel,1) == 1,1,'first')
        resptest_response_all(nChannel,1) = 1;
    else
        resptest_response_all(nChannel,1) = 0;
    end
    
    if find(resptest_response_single(:,nChannel,2) == 2,1,'first')
        resptest_response_all(nChannel,2) = 2;
    else
        resptest_response_all(nChannel,2) = 0;
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