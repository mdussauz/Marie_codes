function [modulated_units,modulation_score, ResidualsMean, ResidualsMedian, ResidualsCI95] = ...
    GetReplayModulatedUnits_v2(SessionName, whichunit, whichARrep, whichPRrep, SigTest)
%version 2 
% - this new version includes the option to pick which repeats to be
%selected
% - it also includes significance testing based on CI95

%INPUT: 
% - SessionName = 'O3/O3_20211005_r0_processed.mat'
% - whichunit = which units to include in analysis / leave empty if user 
% wants to include all cells
% - whichARrep and whichPRrep - eg [4 ; 5] % first for AR repeats and 2nd for PR
% repeats - this corresponds to idx in the OpenLoopPSTH (idx1 of PSTH being the CL
% template) - if one is specified the other has to be specified too
% - SigTest = 'CI' or 'ranksum' - what significance test to use


% OUTPUT:
% - modulated_units = 1 if unit is modulated / 0 if unmodulated
% - modulation_score dim = odor x unit x comparison ((1)CL/AR and (2)CL/PR)
% - ResidualsMean and ResidualsCI95 of all cells to allow to plot scatter plot
% across all sessions

%% What significance testing
if ~exist('SigTest', 'var') 
    SigTest = 'CI'; % default is confidence interval 
end 

%% paths
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

if exist('whichunit', 'var') &&  ~isempty(whichunit)%if no specified units, take all
  SingleUnits = SingleUnits(whichunit);
end

  NbUnit = size(SingleUnits,2);
%% sort units by tetrodes and get open loop rasters and PSTHs
foo = cell2mat(arrayfun(@(x) [x.tetrode; x.id], SingleUnits, 'UniformOutput', false))';
[~, MyUnits] = sortrows(foo,[1 2]);

[OpenLoopTraces,OpenLoopTimestamps,OpenLoopPSTH,OpenLoopRaster] = ...
        ProcessOpenLoopTrialsMD(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'whichunits', MyUnits, 'PSTHsmooth', 100, ...
        'plotfigures', 0);

if exist('whichARrep', 'var') &&  ~isempty(whichARrep)
    whichrep = [1 whichARrep whichPRrep];
    OpenLoopPSTH = OpenLoopPSTH(whichrep,:,:);
    num_OL = numel(whichARrep);
    num_PR = numel(whichPRrep); 
    reps_per_condition = [1 num_OL num_PR];
else
    num_OL = numel(find(ReplayTTLs.TrialID<=max(TrialInfo.TrialID))); % number of active replays
    num_PR = numel(ReplayTTLs.TrialID) - num_OL; % number of passive replays
    reps_per_condition = [1 num_OL num_PR]; % repeats for each condition: CL, OL, PR
end

if strcmp(SessionName,'S12/S12_20230727_r0_processed.mat') 
    if ~exist('whichARrep', 'var') 
        reps_per_condition = [1 8 4];
    end 
end 
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

%% 
switch SigTest

    case 'ranksum'
% Testing whether a cell residual is significant across conditions using Rank sum test
% Rank sum test doesn't work well with less than 4 repeats so instead will
% use CI95

OdorPSTHResiduals = {PSTHResiduals{2:4}};
comparisons = [3 1; 5 2]; % 3-1 = OL-OL vs OL-CL and 5-2 = PR-PR vs PR-CL

for odor = 1:3
    for unit = 1:NbUnit
        for x = 1:size(comparisons,1)
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

    case 'CI' % CI95 comparison
OdorPSTHResiduals = {PSTHResiduals{2:4}};
comparisons = [3 1; 5 2]; % 3-1 = OL-OL vs OL-CL and 5-2 = PR-PR vs PR-CL

for odor = 1:4
    for unit = 1:NbUnit
        for x = 1:size(comparisons,1)
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

modulation_score = modulation_score(2:4,:,:); % not keeping entry 1 which is entire PSTH

end

%% Get modulated units
% % I'm going to assume that only if the residual of cross conditions
% % comparison that are superior to the residual of within condition are
% % worth keeping 
% 
modulated_units   = NaN(NbUnit,2);

for unit = 1:NbUnit
    for condition = 1:size(comparisons,1) %1) is AR vs 2) is PR
        if any(squeeze(modulation_score(:,unit,condition))==1, 'all')
            modulated_units(unit,condition) = 1;
        else
            modulated_units(unit,condition) = 0;
        end
    end
end