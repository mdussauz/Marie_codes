% Reference code for my analysis 
% This contains both correlation AND residual analysis 

%% paths
if strcmp(computer, 'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    datapath = '/mnt/data/Processed/Behavior/';
end

% Mid Session OL
% MySession = fullfile(datapath,'S1', 'S1_20230329_r0_r2_o0_o1_o2_processed.mat');
% MySession = fullfile(datapath,'S3', 'S3_20230329_r0_r1_o0_o1_processed.mat');
% MySession = fullfile(datapath,'S6', 'S6_20230713_r0_r1_o0_o1_processed.mat');
% MySession = fullfile(datapath,'S6', 'S6_20230713_r0_r1_o0_o1_processed.mat');
% MySession = fullfile(datapath,'S7', 'S7_20230705_r0_r1_o0_o1_processed.mat');
MySession = fullfile(datapath,'S11', 'S11_20230813_r0_r1_o0_o1_processed');
% MySession = fullfile(datapath,'S12', 'S12_20230809_r0_r1_o0_o1_processed.mat');

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

N = size(SingleUnits,2);
MyUnits = (1:N);

%% sort units by tetrodes and get open loop rasters and PSTHs
foo = cell2mat(arrayfun(@(x) [x.tetrode; x.id], SingleUnits, 'UniformOutput', false))';
[~, MyUnits] = sortrows(foo,[1 2]);

[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,OpenLoopRaster] = ...
        ProcessOpenLoopTrialsMD(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'whichunits', MyUnits, 'PSTHsmooth', 100, ...
        'plotfigures', 0);

num_OL = numel(find(ReplayTTLs.TrialID<=max(TrialInfo.TrialID))); % active replays
num_PR = numel(ReplayTTLs.TrialID) - num_OL; % passive replays
reps_per_condition = [1 num_OL num_PR]; % CL, OL, PR

%% Bypassing to select specific repeats
OpenLoopPSTH = OpenLoopPSTH([1 2 3 4 8 9 10],:,:);
num_PR = 3;
reps_per_condition = [1 3 3];
%% using the whole PSTH or just select timepoints eg. specific odor for comparison across conditions

% first get full timeseries values
[PSTHCorrs{1}, CorrTags] = ReplayCrossCorr(OpenLoopPSTH,reps_per_condition); % pair wise correlation of the single rep PSTHs
[PSTHResiduals{1}, ResidualTags] = ReplayResiduals(OpenLoopPSTH,reps_per_condition); % pair wise correlation of the single rep PSTHs

% Parse continuous traces into trials
Trial = OpenLoopTraces(:,6);
Trial(Trial<0) = 0;
Odor = abs(OpenLoopTraces(:,6));
TrialTS =  horzcat( find(diff(Odor)>0), find(diff(Trial)>0), find(diff(Trial)<0)); % Odor ON, Trial ON, Trial OFF
TrialTS(:,4) = Odor(TrialTS(:,2)); % which odor
TrialTS(:,5) = OpenLoopTraces(TrialTS(:,2),7); % which Target Zone

OdorSequence = OpenLoop.TTLs.OdorValve{1}(:,4);
if numel(OdorSequence)>size(TrialTS,1)
    OdorSequence(1,:) = [];
end

% Split the long trace into odor-specific stretches 
% only keep points from this trial's odorstart to next trial's odor start
for whichodor = 1:3
    whichones = find(OdorSequence==whichodor);
    MyIdx = [];
    for i = 1:numel(whichones)
        idx(1) = TrialTS(whichones(i),1); % odor start
        if whichones(i)<size(TrialTS,1)
            idx(2) = TrialTS(whichones(i)+1,1) - 1; % next trial odor start
        else
            idx(2) = TrialTS(whichones(i),3) + SampleRate; % 1 sec post trial off
        end
        MyIdx = horzcat(MyIdx,idx(1):idx(2));
    end

    [PSTHCorrs{1+whichodor}] = ReplayCrossCorr(OpenLoopPSTH(:,MyIdx,:),reps_per_condition);
    [PSTHResiduals{1+whichodor}] = ReplayResiduals(OpenLoopPSTH(:,MyIdx,:),reps_per_condition);
end

% %% analyzing behavior traces
% % isolate indices that correspond to 500 ms stretches before and after each
% % trial start
% pairs = [ones(10,1) (2:11)'];
% for i = 1:size(pairs,1)
%     % full trace
%     %thispairCorr = corrcoef(OpenLoopTraces(:,1,pairs(i,1)),OpenLoopTraces(:,1,pairs(i,2)));
%     %pairs(i,3) = thispairCorr(2,1);
%     % lever snipped
%     for j = 1:size(TrialTS,1) % every trial
% %         thistrialcorr = corrcoef(...
% %                         OpenLoopTraces(TrialTS(j,2)+[0:1:100],1,pairs(i,1)),...
% %                         OpenLoopTraces(TrialTS(j,2)+[0:1:100],1,pairs(i,2)) );
% %         pairs(i,2+j) = thistrialcorr(2,1);
%         snip = TrialTS(j,2)+ [-200:1:200];
%         thistrialresidual = OpenLoopTraces(snip,1,pairs(i,1)) - ...
%                             OpenLoopTraces(snip,1,pairs(i,2)) ;
%         pairs(i,2+j) = sqrt(mean(thistrialresidual.^2));
%     end
%     pairs(i,2+j+1) = mean(pairs(i,(3:2+j)));
% end

%% getting means and errors across pairs of either residuals or corrs
% correlations
U = unique(CorrTags); % various types of comparisons - Cl-OL, OL-OL, OL-PR etc
for i = 1:4
    for x = 1:length(U)
        CorrsMedian{i}(:,x) = median(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        CorrsMean{i}(:,x)   = mean(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        CorrsSTD{i}(:,x)    = std(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')';
        %ts = tinv([0.025  0.975],length(x)-1);      % T-Score
        nsamps              = numel(find(CorrTags==U(x)));
        ts                  = tinv([0.025  0.975],(nsamps-1));
        CorrsCI95{i}(:,x)   = ts(2)*std(PSTHCorrs{i}(find(CorrTags==U(x)),:),'omitnan')'/sqrt(nsamps);
    end
end

% residuals
U = unique(ResidualTags); % various types of comparisons - Cl-OL, OL-OL, OL-PR etc
for i = 1:4
    for x = 1:length(U)
        ResidualsMedian{i}(:,x)     = median(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsMean{i}(:,x)       = mean(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsSTD{i}(:,x)        = std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        nsamps                      = numel(find(ResidualTags==U(x)));
        ts                          = tinv([0.025  0.975],(nsamps-1));
        ResidualsCI95{i}(:,x)       = ts(2)*std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')'/sqrt(nsamps);
    end
end

%% CI95 comparison
OdorPSTHResiduals = {PSTHResiduals{2:4}};
comparisons = [3 1; 5 2]; % 3-1 = OL-OL vs OL-CL and 5-2 = PR-PR vs PR-CL

for odor = 1:4
    for unit = 1:N
        for x = 1:length(comparisons)
            temp_median_control = ResidualsMedian{odor}(unit,comparisons(x,1));
            temp_CI_control = ResidualsCI95{odor}(unit,comparisons(x,1));

            temp_median_tested = ResidualsMedian{odor}(unit,comparisons(x,2));
            temp_CI_tested = ResidualsCI95{odor}(unit,comparisons(x,2));
            
            if (temp_median_control + temp_CI_control) < (temp_median_tested - temp_CI_tested)
                modulation_score(odor,unit,x)=1;

            elseif (temp_median_control - temp_CI_control) > (temp_median_tested + temp_CI_tested)
                modulation_score(odor,unit,x)=2;

            else
                modulation_score(odor,unit,x)=0;

            end
        end
    end
end

% I'm going to assume that only if the residual of cross conditions
% conparison that are superior to the residual of withing condition are
% worth keeping 

modulated_units   = NaN(N,2);

for unit = 1:N
    for condition = 1:2
        if any(squeeze(modulation_score(2:4,unit,condition))==1, 'all')
            modulated_units(unit,condition) = 1;
        else
            modulated_units(unit,condition) = 0;
        end
    end
end
%% Plotting the outcomes
% Reorder Units by decreasing OL-OL correlation
[~,SortedbyCorr] = sort(CorrsMedian{1}(:,3),'descend');
SortedbyCorr = circshift(SortedbyCorr,-1); % first value was NaN'

whichtype = 4;
MedianCorrs = CorrsMedian{whichtype};
STDCorrs = CorrsSTD{whichtype};
MedianResiduals = ResidualsMedian{whichtype};
STDResiduals = ResidualsSTD{whichtype};

figure;
subplot(4,1,1); % correlation CL-AR (black) vs AR-AR (pink) 
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianCorrs(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,3) + STDCorrs(SortedbyCorr,3))'; (MedianCorrs(SortedbyCorr,3) - STDCorrs(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);

%CL-OL
bar(xpts, MedianCorrs(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,1) + STDCorrs(SortedbyCorr,1))'; (MedianCorrs(SortedbyCorr,1) - STDCorrs(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);
title('Correlation AR-AR(pink) vs CL-AR (black)')

subplot(4,1,2); % correlation CL-PR (black) vs PR-PR (green) 
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianCorrs(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,5) + STDCorrs(SortedbyCorr,5))'; (MedianCorrs(SortedbyCorr,5) - STDCorrs(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianCorrs(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,2) + STDCorrs(SortedbyCorr,2))'; (MedianCorrs(SortedbyCorr,2) - STDCorrs(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);
title('Correlation PR-PR(green) vs CL-PR (black)')

subplot(4,1,3); % residual CL-AR (black) vs AR-AR (pink) 
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianResiduals(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,3) + STDResiduals(SortedbyCorr,3))'; (MedianResiduals(SortedbyCorr,3) - STDResiduals(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
%CL-OL
bar(xpts, MedianResiduals(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,1) + STDResiduals(SortedbyCorr,1))'; (MedianResiduals(SortedbyCorr,1) - STDResiduals(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);
title('Residual AR-AR(pink) vs CL-AR (black)')

subplot(4,1,4); % residual CL-PR (black) vs PR-PR (pink) 
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianResiduals(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,5) + STDResiduals(SortedbyCorr,5))'; (MedianResiduals(SortedbyCorr,5) - STDResiduals(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianResiduals(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,2) + STDResiduals(SortedbyCorr,2))'; (MedianResiduals(SortedbyCorr,2) - STDResiduals(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);
title('Residual PR-PR(green) vs CL-PR (black)')

%% scatter plot of self vs. across condition residuals
figure, 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianCorrs = []; STDCorrs = []; MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4
    %whichtype = 4;
    MedianCorrs = [MedianCorrs;  CorrsMedian{whichtype}];
    STDCorrs = [STDCorrs; CorrsSTD{whichtype}];

    MedianResiduals = [MedianResiduals;  ResidualsMedian{whichtype}];
    STDResiduals = [STDResiduals; ResidualsCI95{whichtype}];
end

    subplot(1,2,1)
    [~,sortorder] = sort(MedianResiduals(:,5));
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) + STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) - STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5),'k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,2),'.r');
    title('CL-PR')
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);
    subplot(1,2,2)
    [~,sortorder] = sort(MedianResiduals(:,3)+STDResiduals(:,3));
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) + STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) - STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3),'k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,1),'.r');
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);
    title('OL-PR')

%%
figure;
N = length(MedianCorrs);
axes('ColorOrder',brewermap(N,'Spectral'),'NextPlot','replacechildren')

%
% plot within OL corrs on x axis
x1 = MedianCorrs(:,3) - [STDCorrs(:,3) -STDCorrs(:,3)];
y1 = [MedianCorrs(:,1) MedianCorrs(:,1)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,3) MedianCorrs(:,3)];
y2 = MedianCorrs(:,1) - [STDCorrs(:,1) -STDCorrs(:,1)];
subplot(1,3,1)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('r'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('r'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square


% plot within PR corrs on x axis
x1 = MedianCorrs(:,5) - [STDCorrs(:,5) -STDCorrs(:,5)];
y1 = [MedianCorrs(:,2) MedianCorrs(:,2)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,5) MedianCorrs(:,5)];
y2 = MedianCorrs(:,2) - [STDCorrs(:,2) -STDCorrs(:,2)];
subplot(1,3,2)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('t'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('t'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square

% plot within OL corrs on x axis
x1 = MedianCorrs(:,3) - [STDCorrs(:,3) -STDCorrs(:,3)];
y1 = [MedianCorrs(:,4) MedianCorrs(:,4)];
% plot OL-CL corr on Y axis
x2 = [MedianCorrs(:,3) MedianCorrs(:,3)];
y2 = MedianCorrs(:,4) - [STDCorrs(:,4) -STDCorrs(:,4)];
subplot(1,3,3)
line(x1',y1','Linewidth',1); %,'color',Plot_Colors('p'))
hold on
line(x2',y2','Linewidth',1); %,'color',Plot_Colors('p'))
line([-0.2 1], [-0.2 1], 'color', 'k', 'LineStyle',':');
set(gca,'XLim',[-0.2 1], 'YLim', [-0.2 1]);
axis square


%%
% Plot specific Units
%PlotUnits = [58 35 34 55 21];
PlotUnits = [11 89 20 81 45 8 92 70 57 100 86 42];
ProcessOpenLoopTrialsMD(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'plotephys', 1, 'UnitsPerFig', 5, 'whichunits', PlotUnits, 'PlotOpenLoop', 1);


%% bar plots for selected units
PlotUnits = [21 55 34 44 8 39 10 28];
MedianCorrs = []; STDCorrs = []; MedianResiduals =[]; STDResiduals =[];
for k = 1:length(PlotUnits)
    for whichtype = 1:4
        MedianCorrs = [MedianCorrs; CorrsMedian{whichtype}(PlotUnits(k),:)];
        STDCorrs = [STDCorrs; CorrsSTD{whichtype}(PlotUnits(k),:)];
        MedianResiduals = [MedianResiduals; ResidualsMean{whichtype}(PlotUnits(k),:)];
        STDResiduals = [STDResiduals; ResidualsSTD{whichtype}(PlotUnits(k),:)];
    end
end

N = 1:size(MedianResiduals,1);
SortedbyCorr = 1:size(MedianResiduals,1);

figure;
subplot(4,1,1);
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianCorrs(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,3) + STDCorrs(SortedbyCorr,3))'; (MedianCorrs(SortedbyCorr,3) - STDCorrs(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
%CL-OL
bar(xpts, MedianCorrs(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,1) + STDCorrs(SortedbyCorr,1))'; (MedianCorrs(SortedbyCorr,1) - STDCorrs(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);

subplot(4,1,2);
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianCorrs(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,5) + STDCorrs(SortedbyCorr,5))'; (MedianCorrs(SortedbyCorr,5) - STDCorrs(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianCorrs(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianCorrs(SortedbyCorr,2) + STDCorrs(SortedbyCorr,2))'; (MedianCorrs(SortedbyCorr,2) - STDCorrs(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
%xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);

subplot(4,1,3);
xpts = (1:4:4*N);
% OL-OL
bar(xpts, MedianResiduals(SortedbyCorr,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,3) + STDResiduals(SortedbyCorr,3))'; (MedianResiduals(SortedbyCorr,3) - STDResiduals(SortedbyCorr,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
%CL-OL
bar(xpts, MedianResiduals(SortedbyCorr,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,1) + STDResiduals(SortedbyCorr,1))'; (MedianResiduals(SortedbyCorr,1) - STDResiduals(SortedbyCorr,1))'], ...
    'color', [0 0 0],'Linewidth', 2);

subplot(4,1,4);
% PR-PR
xpts = (2:4:4*N);
bar(xpts, MedianResiduals(SortedbyCorr,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,5) + STDResiduals(SortedbyCorr,5))'; (MedianResiduals(SortedbyCorr,5) - STDResiduals(SortedbyCorr,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
% CL-PR
bar(xpts, MedianResiduals(SortedbyCorr,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyCorr,2) + STDResiduals(SortedbyCorr,2))'; (MedianResiduals(SortedbyCorr,2) - STDResiduals(SortedbyCorr,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
%xticklabels(num2str(MyUnits(SortedbyCorr)));
xtickangle(gca,90);