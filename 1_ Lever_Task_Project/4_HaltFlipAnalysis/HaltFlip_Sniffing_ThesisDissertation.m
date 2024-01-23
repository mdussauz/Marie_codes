%HaltFlip_Sniffing_ThesisDissertation
% sniff detection might be slightly under estimating the rate overall 
% needs to be fixed

Sessions = {'O3/O3_20210927_r0_processed.mat',...
'O8/O8_20220704_r0_processed.mat',... 
'O9/O9_20220702_r1_processed.mat',... 
'S1/S1_20230327_r0_processed.mat',... 
'S6/S6_20230710_r0_processed.mat',... 
'S7/S7_20230608_r0_processed.mat',... 
'S11/S11_20230801_r0_processed.mat',... 
'S12/S12_20230731_r0_processed.mat'}; 

OdorForEachMouse = [1 1 1 1 3 2 2 2];
 for mouse = 1: length(Sessions)
     SessionPath = Sessions{mouse};
     whichOdor = OdorForEachMouse(mouse);
     MouseName = fileparts(Sessions{mouse});

%% DataExtraction
if strcmp(computer, 'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
% load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
%                 'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
%                 'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate] = LoadProcessedDataSession(MySession);

%% sniff trace
%sniffs = TracesOut(:,3);
fband = [0.5 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients

old_data = TracesOut(:,3); % for checking purpose
TracesOut(isinf(TracesOut(:,3)),3) = NaN; % some sniff signals have Inf
TracesOut (:,3) = fillmissing(TracesOut (:,3),'nearest');
TracesOut(:,3) = filtfilt(b,a,TracesOut(:,3)); %apply the filter to x(t)


% find peaks
minpeak = 0.05; %used to be 0.2

[pks,pkloc] = findpeaks(TracesOut(:,3),'MinPeakProminence',minpeak,'MinPeakDistance',SampleRate*0.05);
% find valleys
[vls,vlloc] = findpeaks(-TracesOut(:,3),'MinPeakProminence',minpeak,'MinPeakDistance',SampleRate*0.05);

% Plot peaks and sniff trace for checks
figure;
plot(TracesOut(:,3)); hold on; plot(pkloc,pks,"ro"); hold on; plot(vlloc,-vls,"ko"); hold off

% sanity checks
for i = 1:numel(pks)-1 
    % is there a unique valley
    thispkvalley = intersect(find(vlloc>pkloc(i)),find(vlloc<pkloc(i+1)));
    if numel(thispkvalley)==1
        Loc(i,:) = [pkloc(i) vlloc(thispkvalley)];
    else
        %keyboard;
    end
end

%% Trial Aligned Sniff times?
[AlignedSniffs] =  TrialAlignedSniffTimes(Loc,TracesOut(:,5),SampleRate);

%% Get events
% Get all spikes, all units aligned to trials
[AlignedSpikes, Events] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(TrialInfo.TrialID,2),TrialInfo,MySession);

%% Plot
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

whichUnit =1;
[SRs, BinOffset, whichTZ] = ... 
PlotHaltFlips(whichUnit, whichOdor, AlignedSniffs, Events, TrialInfo, AlignTo);

plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(1,:),'color',Plot_Colors('k'),'Linewidth',2);
hold on
plot((1:size(SRs,2))*0.002+BinOffset/1000,SRs(2,:),'color',Plot_Colors('t'),'Linewidth',2);
set(gca,'XLim',myXlim,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12);
legend('normal trials', 'perturbed trials')
ylabel ('Sniff rate (Hz)')
ylim([0 8])
xlabel('Time (s)')
yticks(0:2:10)
xticks([0 6])

%% Calculate sniff rate during perturbation window
perturbation_window = [abs(BinOffset) abs(BinOffset)+1000]/(0.002*1000);
baseline_rate(mouse) = mean(SRs(1,perturbation_window)); 
perturbation_rate(mouse) = mean(SRs(2,perturbation_window)); 

 end

%% summary plot
figure;
subplot(1,2,1) % for each mouse
y = [baseline_rate; perturbation_rate]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];
ylim([0 10])
yticks(0:2:10)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('normal trials', 'perturbed trials')
ylabel ('Sniff rate (Hz)')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [baseline_rate; perturbation_rate]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 10])
yticks(0:2:10)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('Sniff rate (Hz)')
xticklabels({'normal trials', 'perturbed trials'})