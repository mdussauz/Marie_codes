function [resptest_response_class,channels_perc] = ClassifyOdorResponsiveCells (SessionPath, toplot)
% writing this for closed-loop trials only
%ouptuts are :
% - resptest_response_class = for each unit whether 0 unresponsive, 1
% excited, 2 inhibted, 3 mixed
% - percentage of total cells in each class 
% written by MD

%handles.WhereSession.String = fullfile(SessionPath,WhichSession);

perc_plots = toplot;
%% get the data loaded
%MySession = handles.WhereSession.String;
MySession = SessionPath;
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

%% Getting average responses for each odor for air and odor centered period
window_length = 1550; %bin size is 2ms so multiple this number by 2 to get actual time in ms % used to be 3200 but changed to work across mice

AllMeanCenterAlignedFRs = zeros(3, N, window_length); % odor x unit x window length 
AllMeanTrialAlignedFRs = zeros(3, N, window_length); % odor x unit x window length 
AllCenterAlignedFRs = zeros(N, 3, 12, window_length); % unit x odor x tz x window length 
AllTrialAlignedFRs = zeros(N, 3, 12, window_length); % unit x odor x tz x window length 

for whichodor = 1:3
    for whichunit = 1:N
        AlignTo = 5; % to get aligned to TZ entry 
        [AlignedFRs, AlignedRawSpikeCounts] = AlignedFRtoEvent(handles, whichunit, whichodor, AlignTo); %dim are TZ x time 
        if length(AlignedFRs) < window_length % in case the vector is too short to be stored, add zeros to make it the right length
            AlignedFRs(:,window_length) = 0; 
        end 
        %store each tz responses individually:
        AllCenterAlignedFRs(whichunit,whichodor,1:12,1:window_length) = AlignedFRs(:,1:window_length); 
        %take the mean across TZ:
        MeanCenterAlignedFRs = squeeze(mean(AlignedFRs, 1));
        AllMeanCenterAlignedFRs (whichodor,whichunit,:) = MeanCenterAlignedFRs(1:window_length); %store it in matrix that has dimension odor x units x times 
        clearvars MeanCenterAlignedFRs
        
        
        AlignTo = 1; %Trial ON so I can get FR during ITI
        [TrialAlignedFRs, TrialAlignedRawSpikeCounts] = AlignedFRtoEvent(handles, whichunit, whichodor, AlignTo);
        if length(TrialAlignedFRs) < window_length % in case the vector is too short to be stored, add zeros to make it the right length
            TrialAlignedFRs(:,window_length) = 0; 
        end 
        %store each tz responses individually:
        AllTrialAlignedFRs(whichunit,whichodor,1:12,1:window_length) = TrialAlignedFRs(:,1:window_length); 
        %take the mean across TZ:
        MeanTrialAlignedFRs = mean(TrialAlignedFRs, 1);
        AllMeanTrialAlignedFRs (whichodor,whichunit,:) = MeanTrialAlignedFRs(1:window_length); % repeat store in 'baseline' matrix
        clearvars MeanTrialAlignedFRs
        
    end
end

%% Sanity Check of the 2 matrices
%figure(1); plot(squeeze(AllMeanCenterAlignedFRs(3,6,:))); figure(2); plot(squeeze(AllMeanTrialAlignedFRs(3,6,:)))
%% Defining parameters to analyze responsiveness
% input dataset for Diego has dim time x unit x trial
num_Channel                 = N; 
%bin                         = 25;                                          % Set binning for baseline distribution plots 
% target zone has been entered 100 ms before (at bin 1100 here) 
%reference                   = [600 700; 700 800; 800 900; 900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
%trying reduced windows - gives the same result:
%reference                   = [800 900; 900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
% trying even smaller window - final settings - 600ms (300 bins) around target zone entry on each side: 
reference                   = [900 1000; 1000 1100; 1100 1200; 1200 1300; 1300 1400; 1400 1500];  % It is important to use the same window size for each window you want to evaluate.
% bin 600 is when the trial is ON 
reference_baseline          = [100 200; 200 300; 300 400; 400 500];        % You have to divide your baseline in frame windows equivalent to the ones you want to evaluate.
reference_threshold         = [1 99; 5 95; 10 90];                         % Here you can feed possible thresholds you want to try. To get an idea of which threshold to pick, use 'Responsive_bouton_threshold_test_20221021.m'                                   
threshold                   = 1;                                           % Pick the desired threshold by defining the used row of reference threshold.

%% Frame window averages of baseline periods. 
% Here I divided the baseline in 4 periods with a window frame size equivalent to the ones I use for the response period 
AllTrialAlignedFRs = reshape(AllTrialAlignedFRs, [N,36,window_length]); %unit x (odor x tz) x time
 
baseline_A(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(1,1):reference_baseline(1,2)),3,'omitnan')); % 1st window
baseline_B(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(2,1):reference_baseline(2,2)),3,'omitnan')); % 2nd window
baseline_C(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(3,1):reference_baseline(3,2)),3,'omitnan')); % 3rd window
baseline_D(:,:) = squeeze(mean(AllTrialAlignedFRs(:,:,reference_baseline(4,1):reference_baseline(4,2)),3,'omitnan')); % 4th window

baseline_all = [baseline_A baseline_B baseline_C baseline_D];              % These values are concatenated in one single matrix  
%baseline_all = baseline_all'; % to correct for how I organize things vs Diego - dim are channel x (tz x odor x 4 window)
%% Baseline reference values for each channel.
               
set_threshold                   = reference_threshold(threshold,:); 

baseline_Roi_prctl_low_single   = prctile(baseline_all,set_threshold(1,1),2);    % Set the low threshold  
baseline_Roi_prctl_high_single  = prctile(baseline_all,set_threshold(1,2),2);    % Set the high threshold 

baseline_Roi_mean               = mean(baseline_all,2,'omitnan');

baseline_Roi_median             = median(baseline_all,2,'omitnan');

%% Baseline fistribution plots for each ROI. (Pick 'baseline_plots')
baseline_plots = 0; %to plot or not
bin = 25; %20
switch baseline_plots
    case 1
        f2 = figure(2); 
            for nChannel = 1:num_Channel
                    subplot(ceil(sqrt(num_Channel)),ceil(sqrt(num_Channel)),nChannel);
                        histfit(baseline_all(nChannel,:),bin,'normal');
                            title(num2str(nChannel));
                            box off;
                            set(gca,'TickDir','out');
                            %set(gca,'XLim',xrange1);
%                             if nChannel == 1
%                                 xlabel(xlabel1);
%                                 ylabel(ylabel1);
%                             end     
                            vline(baseline_Roi_mean(nChannel,1),'b--');
                            vline(baseline_Roi_median(nChannel,1),'r--');
                            vline(baseline_Roi_prctl_low_single(nChannel,1),'k--');
                            vline(baseline_Roi_prctl_high_single(nChannel,1),'k--');                                                             
            end
    case 0
        f2 = 0;
end
  

%% Period data extraction for responsiveness test
% Currently I am dividing the after cue period in 5 parts to test if in one
% of those I found any changes. You must pick proper frame ranges by
% modifing the 'reference' array and include or remove more dataset windows
% depending on your data.

window_size = reference(1,2)-reference(1,1);% I shouldn't need the +1 - % here for me would be 100 bins = 200ms 

% % Next I reduce the frame and trial dimensions to one by getting an average value

% what we want is a matrix of dimensions unit x
% tested windows which contains average responses to one odor for that
% specific window 
% AllMeanCenterAlignedFRs (whichodor,whichunit,:)

nCategory = size(reference,1); %how many time windows to test 

resptest_values = NaN(3,num_Channel,nCategory); % odor x unit x category
for nPeriods = 1:size(reference,1)
    %calculate mean response for defined time window: 
    resptest_values(:,:,nPeriods) = squeeze(mean(AllMeanCenterAlignedFRs (:,:,reference(nPeriods,1):reference(nPeriods,2)),3)); 
end 
%% Identification of responsive boutons, both enhanced and suppressed for all odors - SINGLE THRESHOLD VALUES


num_Category = nCategory; % the number of tested windows  
resptest_index_single       = NaN(3, num_Channel,num_Category); % odor x unit x time window

% Identification of boutons showing enhanced or suppressed responses for each period
for odor = 1:3
    for nCategory = 1:num_Category
        for nChannel = 1:num_Channel
            if resptest_values(odor,nChannel,nCategory) > baseline_Roi_prctl_high_single(nChannel,1)
                resptest_index_single(odor,nChannel,nCategory) = 1;
                
            elseif  resptest_values(odor, nChannel,nCategory) < baseline_Roi_prctl_low_single(nChannel,1)
                resptest_index_single(odor,nChannel,nCategory) = 2;
                
            else
                resptest_index_single(odor, nChannel,nCategory) = 0;
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
        
%% Percentage analysis of classified channels
%perc_plots = 1;

switch perc_plots
    case 1
        channels_u      = sum(resptest_response_class == 0); %unresponsive
        channels_e      = sum(resptest_response_class == 1); %enhanced
        channels_s      = sum(resptest_response_class == 2); %suppressed
        channels_m      = sum(resptest_response_class == 3); %mixed
        channels_all    = num_Channel;

        channels_perc(1,1) = channels_u*100/channels_all;
        channels_perc(1,3) = channels_e*100/channels_all;
        channels_perc(1,4) = channels_s*100/channels_all;
        channels_perc(1,5) = channels_m*100/channels_all;
        channels_perc(1,2) = channels_perc(1,3)+channels_perc(1,4)+channels_perc(1,5);
        
        resp_channels      = channels_e+channels_s+channels_m;
        
        resp_channels_perc(1,1) = channels_e*100/resp_channels;
        resp_channels_perc(1,2) = channels_s*100/resp_channels;
        resp_channels_perc(1,3) = channels_m*100/resp_channels;

        f3 = figure();
            subplot(1,2,1);
                bar(channels_perc(1,1:2))
                    box off;
                    set(gca,'TickDir','out');
                    xticklabels({'unresponsive','responsive'});
                    ylabel('Units (%)')
                    ylim([0 100]);
                    title('Percentage of responsive units');
                    box off
            subplot(1,2,2);
                bar(resp_channels_perc(1,1:3))
                    box off;
                    set(gca,'TickDir','out');
                    xticklabels({'enhanced','suppressed','mixed'});
                    ylabel('Responsive units (%)')
                    ylim([0 100]);
                    title('Percentage of classification of responses');
                    box off
    case 0    
        f3 = 0;
end    
%%




%% FUNCTIONS 


function [AlignedFRs, AlignedRawSpikeCounts] = AlignedFRtoEvent(handles,whichUnit, whichodor, AlignTo)

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
% of note "whichTrials" is for closed loop trials only and do not take into
% account perturbation trials
for TZ = 1:12
    thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1)); %spiketimes for all trials of 1 tz type
    Events2Align = Offset(find(whichTrials(:,2)==TZ),1); %time at which event is happening for all trials of tz type
    [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500); %resulting psth is an average of all trials of tz type
    AlignedFRs(TZ,1:numel(myFR)) = myFR;
    AlignedRawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
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
end
