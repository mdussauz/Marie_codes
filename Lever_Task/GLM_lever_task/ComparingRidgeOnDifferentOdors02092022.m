function [spikeTrace,fullV1,c1,cFull1,fullV2,c2,cFull2] = ComparingRidgeOnDifferentOdors02092022()
%% what session
WhichSession = 'O3_20210918_r0_processed';
%WhichSession = 'O3_20211005_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% Load the relevant variables

load(handles.WhereSession.String, 'Traces', 'PassiveReplayTraces', 'TrialInfo', 'TargetZones', ...
    'startoffset', 'errorflags', 'SampleRate', ...
    'TTLs', 'ReplayTTLs', 'TuningTTLs', 'SingleUnits');

handles.SessionLength.String = num2str(10*ceil(TTLs.Trial(end,2)/10));
handles.NumUnits.String = num2str(size(SingleUnits,2));
handles.SampleRate = SampleRate;

% Get full behavior traces 
FirstTrialinAnalysis =2; % we chose to ignore trial 1 as potential misalignment with OpEphys
[TracesOut] = ConcatenateTraces(Traces, FirstTrialinAnalysis:length(TrialInfo.TrialID), SampleRate*startoffset);

% calculate the timestamp difference between Ephys and Behavior
TrialStart_behavior = TrialInfo.SessionTimestamps(1,2);
TrialStart_Ephys = TTLs.Trial(1,2);
% factor to convert all behavior timestamps to match Ephys
TimestampAdjuster = TrialStart_Ephys - TrialStart_behavior;

%% Create New variable combining odor and motor move
% Odor Trials 
TrialTrace = TracesOut.Trial{1}; % not really trial but odor ON/OFF
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
MotorTraceAdjusted = TracesOut.Motor{1} + 125; % arbitrary value to remove negative numbers and zeros
%MotorTraceAdjusted = MotorTraceAdjusted.*(Odor1+Odor2+Odor3); % removing motor info when odor is off
% --
Odor1Location = MotorTraceAdjusted.*Odor1;
Odor2Location = MotorTraceAdjusted.*Odor2;
Odor3Location = MotorTraceAdjusted.*Odor3;
AirLocation = MotorTraceAdjusted.*Air;

%% filter the thermistor data
SniffsTrace = TracesOut.Sniffs{1}; % Sniffs
% -- adding extra points to allow for the (known) transient of filter at beginning of signal to
% settle
dataToAdd = 15;
StartBuffer = ones(dataToAdd,1)*SniffsTrace(1);
SniffsTrace = [StartBuffer; SniffsTrace];

SampleRate = 500; % Samples/second
nqf = SampleRate/2; % Nyquist freq.
[b,a] = butter(3,[0.1 30]/nqf,'bandpass');   % Butterworth filter
ThermistorFiltered = filter(b,a,SniffsTrace);  % filtez

% cutting extra datapoints added for filtering
ThermistorFiltered = ThermistorFiltered(dataToAdd+1:end);

%% Create fullR = design matrix 
fullR1 = [Odor1 Odor2 Odor3 MotorTraceAdjusted];
fullR1 = bsxfun(@minus, fullR1, mean(fullR1, 1));

fullR2 = [Odor1Location Odor2Location Odor3Location];
fullR2 = bsxfun(@minus, fullR2, mean(fullR2, 1));

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels1 = {'Odor1' 'Odor2' 'Odor3', 'Motor location'};
regLabels2 = {'Odor1 location' 'Odor2 location' 'Odor3 location'};

%index to reconstruct different response kernels
regIdx1 = [
    ones(1,size(Odor1,2))*find(ismember(regLabels1,'Odor1')) ...
    ones(1,size(Odor2,2))*find(ismember(regLabels1,'Odor2')) ...
    ones(1,size(Odor3,2))*find(ismember(regLabels1,'Odor3'))...
    ones(1,size(MotorTraceAdjusted,2))*find(ismember(regLabels1,'Motor location'))];

regIdx2 = [
    ones(1,size(Odor1Location,2))*find(ismember(regLabels2,'Odor1 location')) ...
    ones(1,size(Odor2Location,2))*find(ismember(regLabels2,'Odor2 location')) ...
    ones(1,size(Odor3Location,2))*find(ismember(regLabels2,'Odor3 location'))];


disp('Regressors and labels computed');

regGroups1 = {'Odor1' 'Odor2' 'Odor3' 'Motor location'};
regGroups1{2,1} = {'Odor1'};
regGroups1{2,2} = {'Odor2'};
regGroups1{2,3} = {'Odor3'};
regGroups1{2,4} = {'Motor location'};

regGroups2 = {'Odor1 location' 'Odor2 location' 'Odor3 location'};
regGroups2{2,1} = {'Odor1 location'};
regGroups2{2,2} = {'Odor2 location'};
regGroups2{2,3} = {'Odor3 location'};


%% PSTH - try

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
    FR(:,i) = MakePSTH(SingleUnits(i).trialstart', 0, PSTHWindow);
    FR(:,i) = FR(:,i)/divideFRby;
end

% PSTH require extra data points to get a round PSTH window to run
% so 1st trimming the data based on the double of lever trace bc PSTH is
% 1000 samples/sec and lever is 500 samples/sec
DatapointToKeep = length(Odor1) *2;
finalFR = FR(1:DatapointToKeep,:);
finalFR = downsample(finalFR,2);

%based on Musall's code
spikeTrace = bsxfun(@minus, finalFR, mean(finalFR)); %make zero-mean

disp('SpikeTrace computed');

%% run full model fit for PSTHs
[~, dimBeta1] = ridgeMML(spikeTrace, fullR1, true); %make model fit
fullFit1 = fullR1 * dimBeta1; %fit data
disp('Full model fit 1 completed');

[~, dimBeta2] = ridgeMML(spikeTrace, fullR2, true); %make model fit
fullFit2 = fullR2 * dimBeta2; %fit data
disp('Full model fit 2 completed');

%% cross-validated models - create for different regressors
ridgeFolds = 10;
fullV1 =  crossValModel1(regLabels1);  %cross-validated full model.
allV1 = NaN(size(fullV1,1), size(fullV1,2), 2, size(regGroups1,2), 'single');
for iRegs1 = 1 : size(regGroups1,2)
    allV1(:, :, 1, iRegs1) = crossValModel1(regGroups1{2,iRegs1}(:)'); %regressor-alone model
    allV1(:, :, 2, iRegs1) = crossValModel1(regLabels1(~ismember(regLabels1,regGroups1{2,iRegs1}(:)'))); %regressor-exclude model
    fprintf('CrossVal 1 for %s complete (%g/%g)\n', regGroups1{1,iRegs1}, iRegs1, size(regGroups1,2));
end
disp('Cross-validation completed');

fullV2 =  crossValModel2(regLabels2);  %cross-validated full model.
allV2 = NaN(size(fullV2,1), size(fullV2,2), 2, size(regGroups2,2), 'single');
for iRegs2 = 1 : size(regGroups2,2)
    allV2(:, :, 1, iRegs2) = crossValModel2(regGroups2{2,iRegs2}(:)'); %regressor-alone model
    allV2(:, :, 2, iRegs2) = crossValModel2(regLabels2(~ismember(regLabels2,regGroups2{2,iRegs2}(:)'))); %regressor-exclude model
    fprintf('CrossVal 2 for %s complete (%g/%g)\n', regGroups2{1,iRegs2}, iRegs2, size(regGroups2,2));
end
disp('Cross-validation completed');
%% plotting
for y1 = 1:size(spikeTrace,2)
    cFull1(y1,1) = corr2(spikeTrace(:,y1), fullV1(:,y1))^2;
    for z1 = 1:size(allV1,4)
        cRegs1(y1,1,z1) = corr2(spikeTrace(:,y1), allV1(:,y1,1,z1))^2;
        cRegs1(y1,2,z1) = corr2(spikeTrace(:,y1), allV1(:,y1,2,z1))^2;
    end
end

for y2 = 1:size(spikeTrace,2)
    cFull2(y2,1) = corr2(spikeTrace(:,y2), fullV2(:,y2))^2;
    for z2 = 1:size(allV2,4)
        cRegs2(y2,1,z2) = corr2(spikeTrace(:,y2), allV2(:,y2,1,z2))^2;
        cRegs2(y2,2,z2) = corr2(spikeTrace(:,y2), allV2(:,y2,2,z2))^2;
    end
end

Odor1Reg1 = ismember(regGroups1(1,:), 'Odor1');
Odor2Reg1 = ismember(regGroups1(1,:), 'Odor2');
Odor3Reg1 = ismember(regGroups1(1,:), 'Odor3');
MotorReg1 = ismember(regGroups1(1,:), 'Motor location');

Odor1Reg2 = ismember(regGroups2(1,:), 'Odor1 location');
Odor2Reg2 = ismember(regGroups2(1,:), 'Odor2 location');
Odor3Reg2 = ismember(regGroups2(1,:), 'Odor3 location');


figure(1)
subplot(1,2,1)
[~, c1] = sort(cFull1,'descend');
plot(cFull1(c1), 'linewidth', 2, 'color', 'k'); hold on;
plot(cRegs1(c1, 1, Odor1Reg1), 'linewidth', 2);
plot(cRegs1(c1, 1, Odor2Reg1), 'linewidth', 2);
plot(cRegs1(c1, 1, Odor3Reg1), 'linewidth', 2);
plot(cRegs1(c1, 1, MotorReg1), 'linewidth', 2);

xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
legend('full model','Odor1', 'Odor2', 'Odor3', 'Motor location')

subplot(1,2,2)
[~, c2] = sort(cFull2,'descend');
plot(cFull2(c2), 'linewidth', 2, 'color', 'k'); hold on;
plot(cRegs2(c2, 1, Odor1Reg2), 'linewidth', 2);
plot(cRegs2(c2, 1, Odor2Reg2), 'linewidth', 2);
plot(cRegs2(c2, 1, Odor3Reg2), 'linewidth', 2);

xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
legend('full model','Odor1 location', 'Odor2 location', 'Odor3 location')
%%
figure(2)
subplot(3,1,1)
plot(TrialTrace(1:20000))
subplot(3,1,2)
plot(ThermistorFiltered(1:20000))
subplot(3,1,3)
plot(spikeTrace(1:20000, 64)); hold on; plot(fullV1(1:20000, 64));hold on; plot(fullV2(1:20000, 64));
legend('Real data','Model 1', 'Model 2')

%% nested function
%% - cheat function for 1st model
function [Vm, cBeta, cLabels] =  crossValModel1(cLabels)
        rng(1) % seed for reproducibility
        Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
        randIdx = randperm(size(spikeTrace,1)); %generate randum number index
        foldCnt = floor(size(spikeTrace,1) / ridgeFolds);
        cIdx = ismember(regIdx1, find(ismember(regLabels1,cLabels))); %get index for task regressors
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            if iFolds == 1
                [cRidge, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR1(dataIdx,cIdx), true); %get beta weights and ridge penalty for task only model
            else
                [~, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR1(dataIdx,cIdx), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
            end
            Vm(~dataIdx, :) = (fullR1(~dataIdx,cIdx) * cBeta); %predict remaining data
            
            if rem(iFolds,ridgeFolds/5) == 0
                fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            end
        end
end

%% - cheat function for 2nd model
function [Vm, cBeta, cLabels] =  crossValModel2(cLabels)
        rng(1) % seed for reproducibility
        Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
        randIdx = randperm(size(spikeTrace,1)); %generate randum number index
        foldCnt = floor(size(spikeTrace,1) / ridgeFolds);
        cIdx = ismember(regIdx2, find(ismember(regLabels2,cLabels))); %get index for task regressors
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            if iFolds == 1
                [cRidge, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR2(dataIdx,cIdx), true); %get beta weights and ridge penalty for task only model
            else
                [~, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR2(dataIdx,cIdx), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
            end
            Vm(~dataIdx, :) = (fullR2(~dataIdx,cIdx) * cBeta); %predict remaining data
            
            if rem(iFolds,ridgeFolds/5) == 0
                fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            end
        end
    end
end