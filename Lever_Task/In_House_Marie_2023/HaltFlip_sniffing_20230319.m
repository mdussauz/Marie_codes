%% for O8
SessionPath = 'O8/O8_20220704_r0_processed.mat';

%% DataExtraction
[Paths] = WhichComputer();
datapath = Paths.Local.Behavior_processed;
MySession = fullfile(datapath,SessionPath);

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
%% Plot closed loop vs active halt 
figure;
AlignTo = 6; 
switch AlignTo
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end

%raster plot - example unit response to active halt and close loop
whichUnit = 21; whichOdor = 1;
subplot(3,2,[1 3]); hold on
[FRs, BinOffset, whichTZ] = ... 
PlotHaltFlipsMD(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

%psth plot - example unit response to active halt and close loop
subplot(3,2,5); hold on
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);

%raster plot - SNIFFING
whichUnit = 1;
subplot(3,2,[2 4]); hold on
[SRs, BinOffset, whichTZ] = ... 
PlotHaltFlipsMD(whichUnit, whichOdor, AlignedSniffs, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

%psth plot - SNIFFING 
subplot(3,2,6); hold on
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);

%% Plot active vs passive halts
figure;
AlignTo = 6; 
switch AlignTo
    case {1,2,6}
        myXlim = [-1.2 6];
    case {3,4}
        myXlim = [-5.2 1];
    case 5
        myXlim = [-2.2 5];
end

%raster plot - example unit response to active halt and close loop
whichUnit = 21; whichOdor = 1;
subplot(3,2,[1 3]); hold on
[FRs, BinOffset, whichTZ] = ... 
PlotHaltFlipsMD(whichUnit, whichOdor, AlignedSpikes, Events, TrialInfo, AlignTo);
set(gca,'XLim',myXlim);

%psth plot - example unit response to active halt and close loop
subplot(3,2,5); hold on
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(FRs,2))*0.002+BinOffset/1000,FRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);

%raster plot - SNIFFING
whichUnit = 1;
subplot(3,2,[2 4]); hold on
[SRs, BinOffset, whichTZ] = ... 
PlotHaltFlipsActivevsPassiveMD(whichUnit, whichOdor, ReplayAlignedSniffs, Events, ReplayInfo, AlignTo, AlignedSpikes, Events, TrialInfo);
%PlotHaltFlipsActivevsPassiveMD(whichUnit, whichOdor, ReplayAlignedSniffs, ReplayEvents, ReplayInfo, AlignTo, AlignedSpikes, Events, TrialInfo);
set(gca,'XLim',myXlim);

%psth plot - SNIFFING 
subplot(3,2,6); hold on
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim);