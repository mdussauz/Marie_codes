%Residual_replay_analysis_script

% Residuals contain the following RMSE of the following comparison for each
% unit:
% 1) CL-OL
% 2) CL-PR
% 3) OL-OL
% 4) PR-PR
% 5) OL-PR

% !!! Later after getting Residual Tags !!!
% PR-PR = 5
% OL-PR = 4

% NB: OL being AR

%% paths


%AON
%classic open loop
%SessionName = fullfile('O3','O3_20211005_r0_processed.mat');
%SessionName = fullfile('O8','O8_20220702_r0_processed.mat');
%SessionName = fullfile('O9','O9_20220630_r0_processed.mat');
%SessionName = fullfile('S1','S1_20230314_r0_processed.mat');
%SessionName = fullfile('S3','S3_20230321_r0_processed.mat');
%SessionName = fullfile('S6','S6_20230727_r0_processed.mat'); % bug
%SessionName = fullfile('S7','S7_20230707_r0_processed.mat');
%SessionName = fullfile('S11','S11_20230812_r0_processed.mat');
%SessionName = fullfile('S12','S12_20230727_r0_processed.mat');

%free open loop 
%SessionName = fullfile('S6','S6_20230718_r0_processed.mat');
%SessionName = fullfile('S7','S7_20230622_r0_processed.mat');
SessionName = fullfile('S11','S11_20230805_r0_processed.mat');


%APC
%SessionName = fullfile('Q4','Q4_20221112_r0_processed.mat'); %currently not sorted 
%SessionName = fullfile('Q9','Q9_20221119_r0_processed.mat');%currently not sorted

if strcmp(computer, 'MACI64')
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    ProcessedBehaviorPath = '/mnt/data/Processed/Behavior/';
end

MySession = fullfile(ProcessedBehaviorPath,SessionName);

%% get the processed data loaded and assemble open loop traces
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

Nb_unit = size(SingleUnits,2);

%% sort units by tetrodes and get open loop rasters and PSTHs
foo = cell2mat(arrayfun(@(x) [x.tetrode; x.id], SingleUnits, 'UniformOutput', false))';
[~, MyUnits] = sortrows(foo,[1 2]);

[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,OpenLoopRaster] = ...
        ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'whichunits', MyUnits, 'PSTHsmooth', 100, ...
        'plotfigures', 0);

num_OL = numel(find(ReplayTTLs.TrialID<=TrialInfo.TrialID(end))); % number of active replays
num_PR = numel(ReplayTTLs.TrialID) - num_OL; % number of passive replays
reps_per_condition = [1 num_OL num_PR]; % repeats for each condition: CL, OL, PR

%% using the whole PSTH or just select timepoints eg. specific odor for comparison across conditions
% PSTHResiduals = [{FullPSTH} {PSTHodor1} {PSTHodor2} {PSTHodor3}]

% first get full timeseries values
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

    [PSTHResiduals{1+whichodor}] = ReplayResiduals(OpenLoopPSTH(:,MyIdx,:),reps_per_condition);
end

%% getting means and errors across pairs of residuals 
% Also reorganizing Residuals pairs as: 
% 1) CL-OL
% 2) CL-PR 
% 3) OL-OL
% 4) OL-PR 
% 5) PR-PR 

U = unique(ResidualTags); % various types of comparisons - Cl-OL, OL-OL, OL-PR etc
SortedPSTHResiduals = [];
for i = 1:4 %entire PSTH, odor 1, 2 and 3
    for x = 1:length(U)
        ResidualsMedian{i}(:,x)     = median(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsMean{i}(:,x)       = mean(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        ResidualsSTD{i}(:,x)        = std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')';
        nsamps                      = numel(find(ResidualTags==U(x)));
        ts                          = tinv([0.025  0.975],(nsamps-1));
        ResidualsCI95{i}(:,x)       = ts(2)*std(PSTHResiduals{i}(find(ResidualTags==U(x)),:),'omitnan')'/sqrt(nsamps);
    end
end

%% Testing whether a cell residual is significant across conditions

OdorPSTHResiduals = {PSTHResiduals{2:4}};
comparisons = [3 1; 5 2]; % 3-1 = OL-OL vs OL-CL and 5-2 = PR-PR vs PR-CL

for odor = 1:3
    for unit = 1:Nb_unit
        for x = 1:length(comparisons)
            % auROC and p-value for ranksum test
            ControlVar = PSTHResiduals{odor+1}(find(ResidualTags==U(comparisons(x,1))),unit); %(:,comparisons(x,1));
            TestedVar = PSTHResiduals{odor+1}(find(ResidualTags==U(comparisons(x,2))),unit); %(unit,:,comparisons(x,2));
            [auROC(odor,unit,x), AURp(odor,unit,x)] = RankSumROC(ControlVar,TestedVar);

            % defining if positively or negatively modulated
            if auROC(odor,unit,x)>.5 && AURp(odor,unit,x)<.05
                modulation_score(odor,unit,x)=1;
            elseif auROC(odor,unit,x)<.5 & AURp(odor,unit,x)<.05
                modulation_score(odor,unit,x)=2; % though here since RMSE doesn't give direction of modulation
            else
                modulation_score(odor,unit,x)=0;
            end
        end
    end
end

%% PLOTTING %%
%% Plotting the outcomes for comparison of PSTHs for odor 3 across conditions
% Reorder Units by decreasing OL-OL residuals
whichtype = 4; % = odor 3 comparison 

MedianResiduals = ResidualsMean{whichtype};
STDResiduals = ResidualsSTD{whichtype};

[~,SortedbyRes] = sort(MedianResiduals(:,3),'descend');
SortedbyRes = circshift(SortedbyRes,-1); % first value was NaN'


figure;
subplot(2,1,1); % residual CL-OL (black) vs OL-OL (pink) 
xpts = (1:4:4*Nb_unit);
% OL-OL
bar(xpts, MedianResiduals(SortedbyRes,3),'Facecolor',[1 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyRes,3) + STDResiduals(SortedbyRes,3))'; (MedianResiduals(SortedbyRes,3) - STDResiduals(SortedbyRes,3))'], ...
    'color', [1 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyRes)));
xtickangle(gca,90);
%CL-OL
bar(xpts, MedianResiduals(SortedbyRes,1),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyRes,1) + STDResiduals(SortedbyRes,1))'; (MedianResiduals(SortedbyRes,1) - STDResiduals(SortedbyRes,1))'], ...
    'color', [0 0 0],'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyRes)));
xtickangle(gca,90);


subplot(2,1,2); % residual CL-PR (black) vs PR-PR (pink) 
% PR-PR
xpts = (2:4:4*Nb_unit);
bar(xpts, MedianResiduals(SortedbyRes,5),'Facecolor',[0.4 0.6 0.6],'BarWidth',0.6,'LineStyle','none');
hold on
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyRes,5) + STDResiduals(SortedbyRes,5))'; (MedianResiduals(SortedbyRes,5) - STDResiduals(SortedbyRes,5))'], ...
    'color', [0.1 0.4 0.2], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyRes)));
xtickangle(gca,90);
% CL-PR
bar(xpts, MedianResiduals(SortedbyRes,2),'Edgecolor','k','Facecolor','none','BarWidth',0.75,'Linewidth', 2);
line(repmat(xpts,2,1), ...
    [(MedianResiduals(SortedbyRes,2) + STDResiduals(SortedbyRes,2))'; (MedianResiduals(SortedbyRes,2) - STDResiduals(SortedbyRes,2))'], ...
    'color', [0 0 0], 'Linewidth', 2);
set(gca, 'XTick', xpts);
xticklabels(num2str(MyUnits(SortedbyRes)));
xtickangle(gca,90);

%% scatter plot of self vs. across condition residuals

figure;
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
    MedianResiduals = [MedianResiduals; ResidualsMean{whichtype}];
    STDResiduals = [STDResiduals; ResidualsCI95{whichtype}];
end

    subplot(1,2,1) %CL-PR with respect to PR-PR
    [~,sortorder] = sort(MedianResiduals(:,5));
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) + STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) - STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5),'k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,2),'or');
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);
    
    subplot(1,2,2) %CL-OL with respect to OL-OL
    [~,sortorder] = sort(MedianResiduals(:,3)+STDResiduals(:,3));
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) + STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) - STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3),'k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,1),'or');
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);