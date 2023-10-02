%% for O8 
SessionPath = 'O8/O8_20220704_r0_processed.mat';
HaltUnits = [21 23 26 31 33 45 50 51 52];
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6):
AlignTo = [2 2 2]; 
FRmax = [0 30];

%% DataExtraction
[Paths] = WhichComputer();
datapath = Paths.Local.Behavior_processed;
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
            
[~, ~, TrialInfo] = LoadProcessedDataSession(MySession);  % to get target zone entry time

%% Passive replays
PassiveTracesOut = []; StartStopIdx = 0;
% make sure open loop template trials have been identified
if ~isempty(PassiveReplayTraces) 
    if ~exist('OpenLoop','var')
        OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);
    end
end
%LoadProcessedPassiveSession;
if ~exist('OpenLoop','var')
    OpenLoop = [];
end

%% Trial Aligned spikeTimes
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

[ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,TrialInfo,Events,OpenLoop);

%% Creating a  new field to facilitate use of UnitPlotter function on Halt Replay 
ReplayInfo.Perturbation(find(ReplayInfo.Perturbed(:,1)==0),1) = {[]};
ReplayInfo.Perturbation(find(ReplayInfo.Perturbed(:,1)==1),1) = {'Halt-Flip-Template'};
%%
whichOdor = 1; 
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
                    find(TrialInfo.Odor==whichOdor));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end

%% Initializing 
FR_tz_cl = NaN(3,5000);
FR_tz_hfa = NaN(3,5000);
FR_tz_hfp = NaN(3,5000);
%% Plotting example active halt responses against closed loop responses
myXlim = [-1.2,5];
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    tz = 1;
    % control trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        
        % first plot control trials
        % Spikes
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, 2, trialsdone, trialboxcolor);
        FR_tz_cl(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
    end
    FR_cl(x,:) = mean(FR_tz_cl,1);
    tz = 1;
    % halt trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(-LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, 6, trialsdone, trialboxcolor);
        FR_tz_hfa(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;
        
        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    FR_hfa(x,:) = mean(FR_tz_hfa,1);
end

%% Plotting example active halt responses against passive halt responses
myXlim = [-1.2,5];
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    tz = 1;
    
    % halt trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(-LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, 6, trialsdone, trialboxcolor);
        
        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    

    % replayed halt trials 
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, ReplayAlignedSpikes, ReplayEvents, ReplayInfo, 6, trialsdone, trialboxcolor);
        FR_tz_hfp(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;

        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    FR_hfp(x,:) = mean(FR_tz_hfp,1);
end

%% MEAN FR 
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichplots = x;
    subplot(1,nCols,whichplots); 
    hold on 
    plot((1:size(FR_cl,2))*0.002+BinOffset/1000,FR_cl(x,:),'color',GetLocationColor(LocationTrialStart(1)),'Linewidth',2);
    plot((1:size(FR_hfa,2))*0.002+BinOffset/1000,FR_hfa(x,:),'color',GetLocationColor(-LocationTrialStart(1)),'Linewidth',2);

    set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
    xticks([0 5])
end 

figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichplots = x;
    subplot(1,nCols,whichplots); 
    hold on
    plot((1:size(FR_hfa,2))*0.002+BinOffset/1000,FR_hfa(x,:),'color',GetLocationColor(-LocationTrialStart(1)),'Linewidth',2);
    plot((1:size(FR_hfp,2))*0.002+BinOffset/1000,FR_hfp(x,:),'color',GetLocationColor(LocationTrialStart(1)),'Linewidth',2);

    set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
    xticks([0 5])
end

%%%
%% GET SNIFF
%% get the processed data loaded

load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

[TracesOut, ColNames, TrialInfo, ~, ~, ...
     ~, ~, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = LoadProcessedDataSession(MySession); %this function gets PassiveTracesOut = [] for some reason
                                                                                                               %TrialInfo of this functions contains target zone entry 
%% sniff trace
%sniffs = TracesOut(:,3);
fband = [0.5 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TracesOut(:,3) = filtfilt(b,a,TracesOut(:,3)); %apply the filter to x(t)

% find peaks
[pks,pkloc] = findpeaks(TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);
% find valleys
[vls,vlloc] = findpeaks(-TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);

% sanity checks
for i = 1:numel(pks)-1 
    % is there a unique valley
    thispkvalley = intersect(find(vlloc>pkloc(i)),find(vlloc<pkloc(i+1)));
    if numel(thispkvalley==1)
        Loc(i,:) = [pkloc(i) vlloc(thispkvalley)];
    else
        keyboard;
    end
end

%% trying to process replay sniff traces
for replaytrials = 1:length(PassiveReplayTraces.Sniffs)
    PassiveReplayTraces.Sniffs{1,replaytrials} = filtfilt(b,a,PassiveReplayTraces.Sniffs{1,replaytrials});

    % find peaks
    [pks,pkloc] = findpeaks(PassiveReplayTraces.Sniffs{1,replaytrials},'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);
    % find valleys
    [vls,vlloc] = findpeaks(-PassiveReplayTraces.Sniffs{1,replaytrials},'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);

    % sanity checks
    for i = 1:numel(pks)-1
        % is there a unique valley
        thispkvalley = intersect(find(vlloc>pkloc(i)),find(vlloc<pkloc(i+1)));
        if numel(thispkvalley==1)
            %Loc(i,:) = [pkloc(i) vlloc(thispkvalley)];
            ReplayAlignedSniffs{replaytrials,1}{1,1}(i) = pkloc(i);
            ReplayAlignedSniffs{replaytrials,2}{1,1}(thispkvalley) = vlloc(thispkvalley);
        else
            keyboard;
        end
    end
end

%%
% figure; 
% subplot(1,3,1); 
% histogram(diff(Loc(:,1))); 
% subplot(1,3,2); 
% histogram(diff(Loc(:,2))); 
% subplot(1,3,3); 
% histogram((Loc(:,2)-Loc(:,1))); 

%% Trial Aligned Sniff times?
[AlignedSniffs] =  TrialAlignedSniffTimes(Loc,TracesOut(:,5),SampleRate);
%% PLOT SNIFF

%% Trial Aligned spikeTimes
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

[ReplayAlignedSpikes, ReplayEvents, ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,TrialInfo,Events,OpenLoop);

clearvars AlignedSpikes ReplayAlignedSpikes
AlignedSpikes = AlignedSniffs;
ReplayAlignedSpikes = ReplayAlignedSniffs;
%% Creating a  new field to facilitate use of UnitPlotter function on Halt Replay 
ReplayInfo.Perturbation(find(ReplayInfo.Perturbed(:,1)==0),1) = {[]};
ReplayInfo.Perturbation(find(ReplayInfo.Perturbed(:,1)==1),1) = {'Halt-Flip-Template'};
%%
whichOdor = 1; 
for tz = 1:12
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
                    find(TrialInfo.Odor==whichOdor));
    whichTrials = intersect(whichTrials,...
                    find(TrialInfo.TargetZoneType==tz));
    OdorONPeriod(tz,1) = median(TrialInfo.OdorStart(whichTrials,1));
    foo = cellfun(@(x) median(x(491:500)), Traces.Motor(whichTrials),'UniformOutput', false);
    LocationTrialStart(tz,1) = round(median(cell2mat(foo)));
end

%% Initializing 
FR_tz_cl = NaN(3,5000);
FR_tz_hfa = NaN(3,5000);
FR_tz_hfp = NaN(3,5000);
%% Plotting example active halt responses against closed loop responses
myXlim = [-1.2,5];
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);
whichUnit =1;
for x = 1:nCols
    %whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    tz = 1;
    % control trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        
        % first plot control trials
        % Spikes
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, 2, trialsdone, trialboxcolor);
        FR_tz_cl(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
    end
    FR_cl(x,:) = mean(FR_tz_cl,1);
    tz = 1;
    % halt trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(-LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, 6, trialsdone, trialboxcolor);
        FR_tz_hfa(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;
        
        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    FR_hfa(x,:) = mean(FR_tz_hfa,1);
end

%% Plotting example active halt responses against passive halt responses
myXlim = [-1.2,5];
nCols = numel(HaltUnits);
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    %whichUnit = HaltUnits(x);
    
    trialsdone = 0;
    tz = 1;
    
    % halt trials
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(-LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, AlignedSpikes, Events, TrialInfo, 6, trialsdone, trialboxcolor);
        
        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    

    % replayed halt trials 
    for whichTZ = 1:4:12
        trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, -whichTZ, ReplayAlignedSpikes, ReplayEvents, ReplayInfo, 6, trialsdone, trialboxcolor);
        FR_tz_hfp(tz,1:size(FRs, 2)) = FRs;
        tz = tz+1;

        set(gca, 'XLim', myXlim,'YColor','none', 'TickDir', 'out');
        xticks([0 5])
        
        trialsdone = trialsdone + nTrials;
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);
        
        set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
        xticks([0 5])
    end
    FR_hfp(x,:) = mean(FR_tz_hfp,1);
end

%% MEAN FR 
figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichplots = x;
    subplot(1,nCols,whichplots); 
    hold on 
    plot((1:size(FR_cl,2))*0.002+BinOffset/1000,FR_cl(x,:),'color',GetLocationColor(LocationTrialStart(1)),'Linewidth',2);
    plot((1:size(FR_hfa,2))*0.002+BinOffset/1000,FR_hfa(x,:),'color',GetLocationColor(-LocationTrialStart(1)),'Linewidth',2);

    set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
    xticks([0 5])
end 

figure('Position',[0 0 1000 600]);

for x = 1:nCols
    whichplots = x;
    subplot(1,nCols,whichplots); 
    hold on
    plot((1:size(FR_hfa,2))*0.002+BinOffset/1000,FR_hfa(x,:),'color',GetLocationColor(-LocationTrialStart(1)),'Linewidth',2);
    plot((1:size(FR_hfp,2))*0.002+BinOffset/1000,FR_hfp(x,:),'color',GetLocationColor(LocationTrialStart(1)),'Linewidth',2);

    set(gca, 'XLim', myXlim, 'YLim', FRmax, 'TickDir', 'out');
    xticks([0 5])
end