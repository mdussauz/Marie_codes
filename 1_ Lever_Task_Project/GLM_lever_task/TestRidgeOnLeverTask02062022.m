function [spikeTrace,fullV,c,cFull] = TestRidgeOnLeverTask02062022()
%% what session

WhichSession = 'O3_20210918_r0_processed';
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
MotorTraceAdjusted = TracesOut.Motor{1} + 125; % arbitrary value to remove negative numbers and zeros
% --
Odor1Location = MotorTraceAdjusted.*Odor1;
Odor2Location = MotorTraceAdjusted.*Odor2;
Odor3Location = MotorTraceAdjusted.*Odor3;
AirLocation = MotorTraceAdjusted.*Air;

%% filter the Lever data
LeverTrace = TracesOut.Lever{1}; % Lever trace
% -- adding extra points to allow for the (known) transient of filter at beginning of signal to
% settle
dataToAdd = 15;
StartBuffer = ones(dataToAdd,1)*LeverTrace(1);
LeverTrace = [StartBuffer; LeverTrace];
SampleRate = 500; % Samples/second
fqc = 25; %passband frequency of the filter in hertz
FilteredLever = lowpass(LeverTrace,fqc,SampleRate);

%% Velocity and acceleration 
LeverVelocity = gradient(LeverTrace);
LeverAcceleration = gradient(LeverVelocity);

FilteredLeverVelocity = gradient(FilteredLever);
FilteredLeverAcceleration = gradient(FilteredLeverVelocity);

% cutting extra datapoints added for filtering
LeverTrace = LeverTrace(dataToAdd+1:end);
FilteredLever = FilteredLever(dataToAdd+1:end);
FilteredLeverVelocity = FilteredLeverVelocity(dataToAdd+1:end);
FilteredLeverAcceleration = FilteredLeverAcceleration(dataToAdd+1:end);

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

%% Rescale the data
ThermistorNormalized = ThermistorFiltered - median(ThermistorFiltered);

%% find points in the thermistor that correspond to the valley
% finding positive peaks
[therm_pks_1,therm_locs_1,therm_w_1,therm_p_1] = findpeaks(ThermistorNormalized, 'MinPeakProminence', 0.01, 'MinPeakDistance', 10);
% finding negative peaks
[therm_pks_2,therm_locs_2,therm_w_2,therm_p_2] = findpeaks(-ThermistorNormalized, 'MinPeakProminence', 0.01, 'MinPeakDistance', 10);

%% combine both positive and negative peaks 
therm_pks_combined = sort([therm_pks_1(:)',therm_pks_2(:)']);
therm_locs_combined = sort([therm_locs_1(:)',therm_locs_2(:)']);

%% Interpolating breathing signal based on peaks of inhalation and exhalation
% -- Adding 1st and final value of Thermistor signal otherwise
% interpolation is wrong at the beginning and at the end
x = [1, therm_locs_combined, length(ThermistorNormalized)];
y = [ThermistorNormalized(1); ThermistorNormalized(therm_locs_combined); ThermistorNormalized(end)];

%% keeping signal between 2 values to remove slow oscillations

%-- defining positive peaks as equal to 2 and negative peaks as equal to 1
y_locs_1 = ones(length(therm_locs_1),1)*2;
y_locs_2 = ones(length(therm_locs_2),1);

%-- getting the index of sorted combined peaks
[B,I] = sort([therm_locs_1(:)',therm_locs_2(:)']);

%-- combining the new 1/2 peaks and sorting them based on index I
y_combined = [y_locs_1(:)', y_locs_2(:)'];
y_combined = y_combined(I);

%-- adding start and end of signal for better interpolation after
%scaling between 1 and 2
scaled_y = normalize(y, 'range', [1 2]);
y_combined = [scaled_y(1); y_combined(:); scaled_y(end)];

%-- interpolating between peaks to get an oscillating signal between peaks
xx = 1:1:length(ThermistorNormalized); 
interpolated_thermistor = pchip(x,y_combined,xx);
interpolated_thermistor = interpolated_thermistor';

%% Create fullR = design matrix 
%fullR = [TrialTrace Odor1Location Odor2Location Odor3Location ThermistorFiltered LeverTrace FilteredLeverVelocity];
%fullR = [TrialTrace Odor1Location Odor2Location Odor3Location interpolated_thermistor FilteredLeverVelocity];
fullR = [TrialTrace Odor1Location Odor2Location Odor3Location interpolated_thermistor];
fullR = bsxfun(@minus, fullR, mean(fullR, 1));

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
%regLabels = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff' 'Lever position' 'Lever velocity'};
%regLabels = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff' 'Lever velocity'};
regLabels = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff'};

%index to reconstruct different response kernels
regIdx = [
    ones(1,size(TrialTrace,2))*find(ismember(regLabels,'Trial On/Off')) ...
    ones(1,size(Odor1Location,2))*find(ismember(regLabels,'Odor1 location')) ...
    ones(1,size(Odor2Location,2))*find(ismember(regLabels,'Odor2 location')) ...
    ones(1,size(Odor3Location,2))*find(ismember(regLabels,'Odor3 location'))...
    ones(1,size(ThermistorFiltered,2))*find(ismember(regLabels,'Sniff'))];
%    ones(1,size(FilteredLeverVelocity,2))*find(ismember(regLabels,'Lever velocity'))
%    ones(1,size(FilteredLever,2))*find(ismember(regLabels,'Lever position')) ...


disp('Regressors and labels computed');

% regressor groups
%regGroups = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff' 'Lever position' 'Lever velocity'}; %group names, second row are individual regressor labels
%regGroups = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff' 'Lever velocity'};
%regGroups = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff' 'Lever velocity'};
regGroups = {'Trial On/Off' 'Odor1 location' 'Odor2 location' 'Odor3 location' 'Sniff'};
regGroups{2,1} = {'Trial On/Off'};
regGroups{2,2} = {'Odor1 location'};
regGroups{2,3} = {'Odor2 location'};
regGroups{2,4} = {'Odor3 location'};
regGroups{2,5} = {'Sniff'};
%regGroups{2,6} = {'Lever velocity'};
%regGroups{2,5} = {'Lever velocity'};
%regGroups{2,5} = {'Lever position'};
%regGroups{2,6} = {'Lever velocity'};

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
DatapointToKeep = length(FilteredLever) *2;
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

TrialReg = ismember(regGroups(1,:), 'Trial On/Off');
Odor1Reg = ismember(regGroups(1,:), 'Odor1 location');
Odor2Reg = ismember(regGroups(1,:), 'Odor2 location');
Odor3Reg = ismember(regGroups(1,:), 'Odor3 location');
SniffReg = ismember(regGroups(1,:), 'Sniff');
%LeverPosReg = ismember(regGroups(1,:), 'Lever position');
%LeverVelReg = ismember(regGroups(1,:), 'Lever velocity');

figure
subplot(1,2,1)
[~, c] = sort(cFull,'descend');
plot(cFull(c), 'linewidth', 2, 'color', 'k'); hold on;
plot(cRegs(c, 1, TrialReg), 'linewidth', 2, 'color', 'r');
plot(cRegs(c, 1, Odor1Reg), 'linewidth', 2, 'color', 'g');
plot(cRegs(c, 1, Odor2Reg), 'linewidth', 2, 'color', 'b');
plot(cRegs(c, 1, Odor3Reg), 'linewidth', 2, 'color', 'c');
plot(cRegs(c, 1, SniffReg), 'linewidth', 2, 'color', 'm');
%plot(cRegs(c, 1, LeverPosReg), 'linewidth', 2, 'color', 'y');
%plot(cRegs(c, 1, LeverVelReg), 'linewidth', 2, 'color', 'w');
xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
%legend('full model','Trial','Odor1', 'Odor2', 'Odor3', 'Sniff', 'LeverPos', 'LeverVel')
%legend('full model','Trial','Odor1', 'Odor2', 'Odor3', 'Sniff', 'LeverVel')
%legend('full model','Trial','Odor1', 'Odor2', 'Odor3', 'LeverVel')
legend('full model','Trial','Odor1', 'Odor2', 'Odor3', 'Sniff')

subplot(1,2,2)
fullMean = [mean(cRegs(:, 1, TrialReg)) mean(cRegs(:, 1, Odor1Reg)) ...
            mean(cRegs(:, 1, Odor2Reg)) mean(cRegs(:, 1, Odor3Reg))...
            mean(cRegs(:, 1, SniffReg))];
           % mean(cRegs(:, 1, LeverVelReg))];  
          %  mean(cRegs(:, 1, LeverPosReg))...
fullError = [sem(cRegs(:, 1, TrialReg)) sem(cRegs(:, 1, Odor1Reg)) ...
             sem(cRegs(:, 1, Odor2Reg)) sem(cRegs(:, 1, Odor3Reg))...
             sem(cRegs(:, 1, SniffReg))];
           %  sem(cRegs(:, 1, LeverVelReg))];
          %   sem(cRegs(:, 1, LeverPosReg))...

errorbar(fullMean,fullError,'k-','linestyle','none','lineWidth',3); hold on
bar(fullMean,'FaceColor','g','EdgeColor','k','BarWidth',0.5,'LineWidth',2);

uniqueMean = [mean(cFull-cRegs(:, 2, TrialReg)) mean(cFull-cRegs(:, 2, Odor1Reg))...
              mean(cFull-cRegs(:, 2, Odor2Reg)) mean(cFull-cRegs(:, 2, Odor3Reg))...
              mean(cFull-cRegs(:, 2, SniffReg))];
              %mean(cFull-cRegs(:, 2, LeverVelReg))];
              %mean(cFull-cRegs(:, 2, LeverPosReg))...
uniqueError = [sem(cFull-cRegs(:, 2, TrialReg)) sem(cFull-cRegs(:, 2, Odor1Reg))...
               sem(cFull-cRegs(:, 2, Odor2Reg)) sem(cFull-cRegs(:, 2, Odor3Reg))...
               sem(cFull-cRegs(:, 2, SniffReg))];
               %sem(cFull-cRegs(:, 2, LeverVelReg))];
               %sem(cFull-cRegs(:, 2, LeverPosReg))...

errorbar(-uniqueMean,uniqueError,'k-','linestyle','none','lineWidth',3); hold on
bar(-uniqueMean,'FaceColor',[25 111 61]/255,'EdgeColor','k','BarWidth',0.5,'LineWidth',2);

ax = gca;
set(ax,'xTick',1:size(fullMean,2))
%set(ax,'xTickLabel',{'Trial','Odor1', 'Odor2', 'Odor3', 'Sniff', 'LeverPos', 'LeverVel'})
%set(ax,'xTickLabel',{'Trial','Odor1', 'Odor2', 'Odor3', 'Sniff', 'LeverVel'})
%set(ax,'xTickLabel',{'Trial','Odor1', 'Odor2', 'Odor3', 'LeverVel'})
set(ax,'xTickLabel',{'Trial','Odor1', 'Odor2', 'Odor3', 'Sniff'})
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