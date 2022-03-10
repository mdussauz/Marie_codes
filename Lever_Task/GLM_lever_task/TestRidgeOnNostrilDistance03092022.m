function [spikeTrace,fullV,c,cFull] = TestRidgeOnNostrilDistance03092022()
%% what session

WhichSession = 'O3_20211005_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% Load relevant variables

load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.SampleRate = SampleRate;

%% Load Template and Open Loop trial info 

x = find(~cellfun(@isempty, TrialInfo.Perturbation(:,1)));
u = unique(TrialInfo.Perturbation(x));
handles.PerturbationList.String = u{1};
for y = 2:size(u,1)
    handles.PerturbationList.String = [handles.PerturbationList.String,'; ',u{y}];
end

% get start and stop TS of the template
templateStart = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'first'));
templateEnd   = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Template'),1,'last'));

% get start and stop of the replays
replays = x(find(strcmp(TrialInfo.Perturbation(x),'OL-Replay')));

%% Get full behavior traces for Baseline Closed Loop 
FirstTrialinAnalysis =2; % we chose to ignore trial 1 as potential misalignment with OpEphys
LastClosedLoopBaseline = replays(1)-1;
[TracesOut] = ConcatenateTraces(Traces, FirstTrialinAnalysis:LastClosedLoopBaseline, SampleRate*startoffset);

%% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

%% Create New variable combining odor and motor move
% Odor Trials 
TrialTrace = TracesOut.Trial{1};
Odor1 = TrialTrace;
Odor1(Odor1~=1) = 0;
Odor2 = TrialTrace;
Odor2(Odor2~=2) = 0;
Odor2(Odor2==2) = 1;
Odor3 = TrialTrace;
Odor3(Odor3~=3) = 0;
Odor3(Odor3==3) = 1;
Air = TrialTrace;
Air(Air~=4) = 0;
Air(Air==4) = 1;
% --
MotorTrace = TracesOut.Motor{1}; % Motor trace - odor movement

% %% filter the Motor data -- seem to make the model fail...
% % -- adding extra points to allow for the (known) transient of filter at beginning of signal to
% % settle
% dataToAdd = 15;
% StartBuffer = ones(dataToAdd,1)*MotorTrace(1);
% MotorTrace = [StartBuffer; MotorTrace];
% SampleRate = 500; % Samples/second
% fqc = 25; %passband frequency of the filter in hertz
% FilteredMotor = lowpass(MotorTrace,fqc,SampleRate);
% MotorTrace = FilteredMotor(dataToAdd+1:end);

%% create lign for each faske nostril and calculate distance
x_length = length(MotorTrace);
rightNostril = ones(x_length,1)*(-5);
leftNostril = ones(x_length,1)*5;

odorDist_rightNostril = abs(MotorTrace - rightNostril);
odorDist_leftNostril = abs(MotorTrace - leftNostril);

%% Creating tine shifted array 
odorDistRightShifted{1} = odorDist_rightNostril;
odorDistLeftShifted{2} = odorDist_leftNostril;

for i = 1:25 
    odorDistRightShifted{i+1} = circshift(odorDist_rightNostril,i,1);
    odorDistLeftShifted{i+1} = circshift(odorDist_leftNostril,i,1);
end 
clear odorDist_rightNostril; clear odorDist_leftNostril;
odorDist_rightNostril = cat(2,odorDistRightShifted{:});
odorDist_leftNostril = cat(2,odorDistLeftShifted{:});

%% Combining with od on/off

odor1Dist_right = odorDist_rightNostril .*Odor1;
odor2Dist_right = odorDist_rightNostril .*Odor2;
odor3Dist_right = odorDist_rightNostril .*Odor3;

odor1Dist_left = odorDist_leftNostril .*Odor1;
odor2Dist_left = odorDist_leftNostril .*Odor2;
odor3Dist_left = odorDist_leftNostril .*Odor3;


%% Create fullR = design matrix 

fullR = [odor1Dist_right odor2Dist_right odor3Dist_right odor1Dist_left odor2Dist_left odor3Dist_left];

fullR = bsxfun(@minus, fullR, mean(fullR, 1));

regLabels = {'Odor1 R' 'Odor2 R' 'Odor3 R' 'Odor1 L' 'Odor2 L' 'Odor3 L'};

%index to reconstruct different response kernels
regIdx = [
    ones(1,size(odor1Dist_right,2))*find(ismember(regLabels,'Odor1 R')) ...
    ones(1,size(odor2Dist_right,2))*find(ismember(regLabels,'Odor2 R')) ...
    ones(1,size(odor3Dist_right,2))*find(ismember(regLabels,'Odor3 R'))...
    ones(1,size(odor1Dist_left,2))*find(ismember(regLabels,'Odor1 L'))...
    ones(1,size(odor2Dist_left,2))*find(ismember(regLabels,'Odor2 L'))...
    ones(1,size(odor3Dist_left,2))*find(ismember(regLabels,'Odor3 L'))...
    ];

disp('Regressors and labels computed');

% regressor groups
regGroups = {'Odor1 R' 'Odor2 R' 'Odor3 R' 'Odor1 L' 'Odor2 L' 'Odor3 L'};;
regGroups{2,1} = {'Odor1 R'};
regGroups{2,2} = {'Odor2 R'};
regGroups{2,3} = {'Odor3 R'};
regGroups{2,4} = {'Odor1 L'};
regGroups{2,5} = {'Odor2 L'};
regGroups{2,6} = {'Odor3 L'};

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
    %FR(:,i) = MakePSTH(SingleUnits(i).trialstart', 0, PSTHWindow, 'kernelsize',25);
    FR(:,i) = MakePSTH(SingleUnits(i).trialstart', 0, PSTHWindow);
    FR(:,i) = FR(:,i)/divideFRby;
    %newFR(:,i) = downsample(FR(:,i),100); %this gives a less smooth PSTH
end

% PSTH require extra data points to get a round PSTH window to run
% so 1st trimming the data based on the double of lever trace bc PSTH is
% 1000 samples/sec and lever is 500 samples/sec
DatapointToKeep = length(fullR) *2;
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

Odor1RReg = ismember(regGroups(1,:), 'Odor1 R');
Odor2RReg = ismember(regGroups(1,:), 'Odor2 R');
Odor3RReg = ismember(regGroups(1,:), 'Odor3 R');
Odor1LReg = ismember(regGroups(1,:), 'Odor1 L');
Odor2LReg = ismember(regGroups(1,:), 'Odor2 L');
Odor3LReg = ismember(regGroups(1,:), 'Odor3 L');


figure(1)
subplot(1,2,1)
[~, c] = sort(cFull,'descend');
plot(cFull(c), 'linewidth', 2, 'color', 'k'); hold on;
plot(cRegs(c, 1, Odor1RReg), 'linewidth', 2);
plot(cRegs(c, 1, Odor2RReg), 'linewidth', 2);
plot(cRegs(c, 1, Odor3RReg), 'linewidth', 2);
plot(cRegs(c, 1, Odor1LReg), 'linewidth', 2);
plot(cRegs(c, 1, Odor2LReg), 'linewidth', 2);
plot(cRegs(c, 1, Odor3LReg), 'linewidth', 2);

xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
legend('full model','Odor1 R', 'Odor2 R', 'Odor3 R', 'Odor1 L', 'Odor2 L', 'Odor3 L')

%%
%%
figure(2)
subplot(3,1,1)
plot(TrialTrace(1:20000), 'color', 'k')
subplot(3,1,2)
whatneuron = c(1); 
plot(spikeTrace(1:20000, whatneuron),'k'); hold on; plot(fullV(1:20000, whatneuron),'r');
legend('Real data','Model')
subplot(3,1,3)
whatneuron = c(2); 
plot(spikeTrace(1:20000, whatneuron),'k'); hold on; plot(fullV(1:20000, whatneuron),'r');
legend('Real data','Model')

%% nested function
function [Vm, cBeta, cLabels] =  crossValModel(cLabels)
        %rng(1) % for reproducibility
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