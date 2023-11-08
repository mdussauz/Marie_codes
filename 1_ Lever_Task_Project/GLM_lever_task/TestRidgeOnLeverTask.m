function [spikeTrace,fullV,c,cFull] = Test()
%% what session
%WhichSession = 'O1_20211021_r0_processed';
WhichSession = 'O3_20210918_r0_processed';
% SessionPath = 'Z:\mdussauz\Smellocator\Processed\Behavior\O3';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% Load the relevant variables

load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));

% Perturbations
if any(~cellfun(@isempty, TrialInfo.Perturbation(:,1)))
    x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
    u = unique(TrialInfo.Perturbation(x));
    handles.PerturbationList.String = u{1};
    for y = 2:size(u,1)
        handles.PerturbationList.String = [handles.PerturbationList.String,'; ',u{y}];
    end
else
    handles.PerturbationList.String = '';
end

handles.SampleRate = SampleRate;

% Get full behavior traces 
FirstTrialinAnalysis =2; % we chose to ignore trial 1 as potential misalignment with OpEphys
%[TracesOut] = ConcatenateTraces(Traces, 1:length(TrialInfo.TrialID), SampleRate*startoffset);
[TracesOut] = ConcatenateTraces(Traces, FirstTrialinAnalysis:length(TrialInfo.TrialID), SampleRate*startoffset);

% create a corresponding timestamp vector
Timestamps = TracesOut.Timestamps{1}';

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

% Lever trace
LeverTrace = TracesOut.Lever{1};

% Motor trace - odor movement
MotorTrace = TracesOut.Motor{1};

% Rewards
rewards_timestamps =  Timestamps(find(diff(TracesOut.Rewards{1}==1)) + 1)' + TimestampAdjuster;
RewardsTrace = TracesOut.Rewards{1};

% Sniffs
SniffsTrace = TracesOut.Sniffs{1};

% Licks
LicksTrace = TracesOut.Licks{1};

% Odor Trials 
TrialTrace = TracesOut.Trial{1};
Odor1 = TrialTrace;
Odor1(Odor1~=1) = 0;
Odor2 = TrialTrace;
Odor2(Odor2~=2) = 0;
Odor3 = TrialTrace;
Odor3(Odor3~=3) = 0;
Air = TrialTrace;
Air(Air~=4) = 0;

%% Create fullR = design matrix 

switch WhichSession
    case 'O3_20210918_r0_processed_ignore'
        % fullR = [LeverTrace MotorTrace RewardsTrace SniffsTrace LicksTrace Odor1 Odor2 Odor3 Air];
        fullR = [LeverTrace MotorTrace RewardsTrace SniffsTrace LicksTrace Odor1 Odor2 Odor3];
        fullR = bsxfun(@minus, fullR, mean(fullR, 1));
        
        % labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
        %regLabels = {'Lever' 'Motor' 'Rewards' 'Sniff' 'Licks' 'Odor1' 'Odor2' 'Odor3' 'Air'};
        regLabels = {'Lever' 'Motor' 'Rewards' 'Sniff' 'Licks' 'Odor1' 'Odor2' 'Odor3'};
        
        %index to reconstruct different response kernels
        regIdx = [
            ones(1,size(LeverTrace,2))*find(ismember(regLabels,'Lever')) ...
            ones(1,size(MotorTrace,2))*find(ismember(regLabels,'Motor')) ...
            ones(1,size(RewardsTrace,2))*find(ismember(regLabels,'Rewards')) ...
            ones(1,size(SniffsTrace,2))*find(ismember(regLabels,'Sniff')) ...
            ones(1,size(LicksTrace,2))*find(ismember(regLabels,'Licks')) ...
            ones(1,size(Odor1,2))*find(ismember(regLabels,'Odor1')) ...
            ones(1,size(Odor2,2))*find(ismember(regLabels,'Odor2')) ...
            ones(1,size(Odor3,2))*find(ismember(regLabels,'Odor3'))];
        %     ones(1,size(Odor3,2))*find(ismember(regLabels,'Odor3'))...
        %     ones(1,size(Air,2))*find(ismember(regLabels,'Air'))];
        
        disp('Regressors and labels computed');
        
        % regressor groups
        regGroups = {'Lever' 'Motor' 'Rewards' 'Sniff' 'Licks' 'Odor'}; %group names, second row are individual regressor labels
        %regGroups = {'Lever' 'Motor' 'Rewards' 'Sniff' 'Licks' 'Odor' 'Air' 'centerSpout'}; %group names, second row are individual regressor labels
        regGroups{2,1} = {'Lever'};
        regGroups{2,2} = {'Motor'};
        regGroups{2,3} = {'Rewards'};
        regGroups{2,4} = {'Sniff'};
        regGroups{2,5} = {'Licks'};
        regGroups{2,6} = {'Odor1' 'Odor2' 'Odor3'};
        % regGroups{2,7} = {'Air'};
        % regGroups{2,8} = {'Odor1' 'Odor2' 'Odor3' 'Air'};
        
    otherwise
        fullR = [LeverTrace MotorTrace Odor1 Odor2 Odor3 SniffsTrace];
        fullR = bsxfun(@minus, fullR, mean(fullR, 1));
        
        % labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
        regLabels = {'Lever' 'Motor' 'Odor1' 'Odor2' 'Odor3' 'Sniff'};
        
        %index to reconstruct different response kernels
        regIdx = [
            ones(1,size(LeverTrace,2))*find(ismember(regLabels,'Lever')) ...
            ones(1,size(MotorTrace,2))*find(ismember(regLabels,'Motor')) ...
            ones(1,size(Odor1,2))*find(ismember(regLabels,'Odor1')) ...
            ones(1,size(Odor2,2))*find(ismember(regLabels,'Odor2')) ...
            ones(1,size(Odor3,2))*find(ismember(regLabels,'Odor3'))...
            ones(1,size(SniffsTrace,2))*find(ismember(regLabels,'Sniff'))];
        
        disp('Regressors and labels computed');
        
        % regressor groups
        regGroups = {'Lever' 'Motor' 'Odor1' 'Odor2' 'Odor3' 'Sniff'}; %group names, second row are individual regressor labels
        regGroups{2,1} = {'Lever'};
        regGroups{2,2} = {'Motor'};
        regGroups{2,3} = {'Odor1'};
        regGroups{2,4} = {'Odor2'};
        regGroups{2,5} = {'Odor3'};
        regGroups{2,6} = {'Sniff'};
end

%% PSTH - try
% SessionLength = str2num(handles.SessionLength.String);
%PSTHWindow = [0 1000*SessionLength];

downsamplefactor = 100;
divideFRby = 20;

% Remove spikes that happened before 2nd trial and after last trial
% And adjust OpEphys timestamps to match behavior data
SessionStart = TTLs.Trial(2,1)-1; % lever trace starts 1 sec before trial start + starting at trial 2
SessionEnd = TTLs.Trial(end,2)+1 ; % very end of lever trace ends 1 sec after trial start % need to check with Priyanka

for myUnit = 1:length(SingleUnits) % for each cluster
    
    allspikes = SingleUnits(myUnit).spikes ; % in seconds % align to behavior
    SingleUnits(myUnit).trialstart = ...
            SingleUnits(myUnit).spikes(find(allspikes>=SessionStart & allspikes<=SessionEnd)) - SessionStart  ;
end

PSTHendWindow = 100*ceil(SessionEnd/100);
PSTHWindow = [0 1000*PSTHendWindow];

for i = 1:1:size(SingleUnits,2)
    %FR(:,i) = MakePSTH(SingleUnits(i).trialstart', 0, PSTHWindow, 'downsample', downsamplefactor);
    FR(:,i) = MakePSTH(SingleUnits(i).trialstart', 0, PSTHWindow, 'kernelsize',100);
    FR(:,i) = FR(:,i)/divideFRby;
    %newFR(:,i) = downsample(FR(:,i),100); %this gives a less smooth PSTH
end

% PSTH require extra data points to get a round PSTH window to run
% so 1st trimming the data based on the double of lever trace bc PSTH is
% 1000 samples/sec and lever is 500 samples/sec
DatapointToKeep = length(LeverTrace) *2;
finalFR = FR(1:DatapointToKeep,:);
finalFR = downsample(finalFR,2);

%based on Musall's code
spikeTrace = bsxfun(@minus, finalFR, mean(finalFR)); %make zero-mean

disp('SpikeTrace computed');

%% run full model fit for PSTHs
[~, dimBeta] = ridgeMML(spikeTrace, fullR, true); %make model fit
fullFit = fullR * dimBeta; %fit data
% stimIdx = regIdx == find(ismember(regLabels,'stim'));
% % vidFit = fullR(:, ~stimIdx) * dimBeta(~stimIdx, :); %fit data
% stimFit = fullR(:, stimIdx) * dimBeta(stimIdx, :); %fit data
% save([fPath 'fullFit.mat'], 'fullFit', 'vidFit', 'stimFit');
% save([fPath 'regData.mat'], 'fullR', 'dimBeta', 'regLabels', 'regIdx');
disp('Full model fit completed');

%% cross-validated models - create for different regressors
ridgeFolds = 10;
fullV =  crossValModel(regLabels);  %cross-validated full model.
allV = NaN(size(fullV,1), size(fullV,2), 2, size(regGroups,2), 'single');
for iRegs = 1 : size(regGroups,2)
    allV(:, :, 1, iRegs) = crossValModel(regGroups{2,iRegs}(:)'); %regressor-alone model
    allV(:, :, 2, iRegs) = crossValModel(regLabels(~ismember(regLabels,regGroups{2,iRegs}(:)'))); %regressor-exclude model
    fprintf('CrossVal for %s complete (%g/%g)\n', regGroups{1,iRegs}, iRegs, size(regGroups,2));
end
disp('Cross-validation completed');

%% plotting
for y = 1:size(spikeTrace,2)
    cFull(y,1) = corr2(spikeTrace(:,y), fullV(:,y))^2;
    for z = 1:size(allV,4)
        cRegs(y,1,z) = corr2(spikeTrace(:,y), allV(:,y,1,z))^2;
        cRegs(y,2,z) = corr2(spikeTrace(:,y), allV(:,y,2,z))^2;
    end
end

% cFull = cat(1,cFull{:});
% cRegs = cat(1,cRegs{:});
% % concatenating here might not be necessary 

Odor1Reg = ismember(regGroups(1,:), 'Odor1');
Odor2Reg = ismember(regGroups(1,:), 'Odor2');
Odor3Reg = ismember(regGroups(1,:), 'Odor3');
LeverReg = ismember(regGroups(1,:), 'Lever');
MotorReg = ismember(regGroups(1,:), 'Motor');
SniffReg = ismember(regGroups(1,:), 'Sniff');
% RewardsReg = ismember(regGroups(1,:), 'Rewards');
% LicksReg = ismember(regGroups(1,:), 'Licks');

figure
subplot(1,2,1)
[~, c] = sort(cFull,'descend');
plot(cFull(c), 'linewidth', 4, 'color', 'k'); hold on;
plot(cRegs(c, 1, Odor1Reg), 'linewidth', 2);
plot(cRegs(c, 1, Odor2Reg), 'linewidth', 2);
plot(cRegs(c, 1, Odor3Reg), 'linewidth', 2);
plot(cRegs(c, 1, LeverReg), 'linewidth', 2);
plot(cRegs(c, 1, MotorReg), 'linewidth', 2);
xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
legend('full model','Odor1','Odor2','Odor3', 'Lever', 'Motor')

subplot(1,2,2)
% fullMean = [mean(cRegs(:, 1, OdorReg)) mean(cRegs(:, 1, LeverReg)) mean(cRegs(:, 1, MotorReg))...
%     mean(cRegs(:, 1, SniffReg)) mean(cRegs(:, 1, RewardsReg)) mean(cRegs(:, 1, LicksReg))];
% fullError = [sem(cRegs(:, 1, OdorReg)) sem(cRegs(:, 1, LeverReg)) sem(cRegs(:, 1, MotorReg))...
%     sem(cRegs(:, 1, SniffReg)) sem(cRegs(:, 1, RewardsReg)) sem(cRegs(:, 1, LicksReg))];
fullMean = [mean(cRegs(:, 1, OdorReg)) mean(cRegs(:, 1, LeverReg)) mean(cRegs(:, 1, MotorReg)) mean(cRegs(:, 1, SniffReg))];
fullError = [sem(cRegs(:, 1, OdorReg)) sem(cRegs(:, 1, LeverReg)) sem(cRegs(:, 1, MotorReg)) sem(cRegs(:, 1, SniffReg))];

errorbar(fullMean,fullError,'k-','linestyle','none','lineWidth',3); hold on
bar(fullMean,'FaceColor','g','EdgeColor','k','BarWidth',0.5,'LineWidth',2);

% uniqueMean = [mean(cFull-cRegs(:, 2, OdorReg)) mean(cFull-cRegs(:, 2, LeverReg)) mean(cFull-cRegs(:, 2, MotorReg))...
%     mean(cFull-cRegs(:, 2, SniffReg)) mean(cFull-cRegs(:, 2, RewardsReg)) mean(cFull-cRegs(:, 2, LicksReg))];
% uniqueError = [sem(cFull-cRegs(:, 2, OdorReg)) sem(cFull-cRegs(:, 2, LeverReg)) sem(cFull-cRegs(:, 2, MotorReg))...
%     sem(cFull-cRegs(:, 2, SniffReg)) sem(cFull-cRegs(:, 2, RewardsReg)) sem(cFull-cRegs(:, 2, LicksReg))];
uniqueMean = [mean(cFull-cRegs(:, 2, OdorReg)) mean(cFull-cRegs(:, 2, LeverReg)) mean(cFull-cRegs(:, 2, MotorReg)) mean(cFull-cRegs(:, 2, SniffReg))];
uniqueError = [sem(cFull-cRegs(:, 2, OdorReg)) sem(cFull-cRegs(:, 2, LeverReg)) sem(cFull-cRegs(:, 2, MotorReg)) sem(cFull-cRegs(:, 2, SniffReg))];

errorbar(-uniqueMean,uniqueError,'k-','linestyle','none','lineWidth',3); hold on
bar(-uniqueMean,'FaceColor',[25 111 61]/255,'EdgeColor','k','BarWidth',0.5,'LineWidth',2);

ax = gca;
set(ax,'xTick',1:size(fullMean,2))
set(ax,'xTickLabel',{'Odor' 'Lever', 'Motor', 'Sniff', 'Rewards', 'Licks'})
set(ax,'XTickLabelRotation',45)
ax.TickLength = [0 0];
ylabel('cross-val. R^2 R^2'); ylim([-0.2 0.2]);
axis square

%% nested function
function [Vm, cBeta, cLabels] =  crossValModel(cLabels)
        rng(1) % for reproducibility
        Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
        randIdx = randperm(size(spikeTrace,1)); %generate randum number index
        foldCnt = floor(size(spikeTrace,1) / ridgeFolds);
        cIdx = ismember(regIdx, find(ismember(regLabels,cLabels))); %get index for task regressors
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            if iFolds == 1
                [cRidge, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true); %get beta weights and ridge penalty for task only model
            else
                [~, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
            end
            Vm(~dataIdx, :) = (fullR(~dataIdx,cIdx) * cBeta); %predict remaining data
            
            if rem(iFolds,ridgeFolds/5) == 0
                fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            end
        end
    end

end