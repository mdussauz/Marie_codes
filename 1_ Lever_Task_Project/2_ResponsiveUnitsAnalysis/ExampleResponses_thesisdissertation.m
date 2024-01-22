% %% for O3
% SessionPath = 'O3/O3_20210927_r0_processed.mat';
% % Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
% AlignTo = [2 2 2]; 
% ChosenUnits = [18 25 53];
% FRmax = [0 30];
 
% for S6
SessionPath = 'S6/S6_20230710_r0_processed.mat';
% Trial ON (1), odor ON (2), Trial OFF (3), Reward (4), Halt Start (6)
AlignTo = [1 1 1 1 1 1];
ChosenUnits = [1 88 32 83 82 30];
whichOdor = 3;
FRmax = [0 30];
LocationTrialStart = -100;
TZ = 1;

% 1/87 = trial one
% 88 = excited by all odors
% 83 = excited by odor 2 and 3
% 38/82 = excited by odor 3
% 24/32 = trial one
% 30 = odor 1


%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
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

%% Plot
nCols = numel(ChosenUnits);
figure('Position',[0 0 1000 600]);
for x = 1:nCols
    whichUnit = ChosenUnits(x);
    
    switch AlignTo(x)
        case {1,2,6}
            myXlim = [-1 5];
        case {3,4}
            myXlim = [-5 1];
        case 5
            myXlim = [-1.2 5];
    end

    trialsdone = 0;
    for whichTZ = TZ
        
        %trialboxcolor = GetLocationColor(LocationTrialStart(whichTZ));
        trialboxcolor = GetLocationColor(LocationTrialStart);
        
        % Spikes
        whichplots = x + [0 nCols];
        subplot(3,nCols,whichplots);
        hold on
        [nTrials, FRs, BinOffset] = ...
            UnitPlotter_v2(whichUnit, whichOdor, whichTZ, AlignedSpikes, Events, TrialInfo, AlignTo(x), trialsdone, trialboxcolor);
        
        set(gca, 'XLim', myXlim, 'YLim', [0 7], 'TickDir', 'out');
        
        % FR
        whichplots = x + 2*nCols;
        subplot(3,nCols,whichplots); 
        hold on
        plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',trialboxcolor,'Linewidth',2);

        set(gca, 'XLim', myXlim, 'YLim', [0 max(FRs)+5], 'TickDir', 'out');
        
        trialsdone = trialsdone + nTrials;
    end    
end