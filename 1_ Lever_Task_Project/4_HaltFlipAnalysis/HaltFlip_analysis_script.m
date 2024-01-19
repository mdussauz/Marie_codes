%HaltFlip_analysis_script

%% Sessions
%AON
%SessionName = 'O3/O3_20210927_r0_processed.mat';
%SessionName = 'O3/O3_20210929_r0_processed.mat';
SessionName = 'O8/O8_20220704_r0_processed.mat';
%SessionName = 'O9/O9_20220702_r1_processed.mat';
%SessionName = 'S1/S1_20230327_r0_processed.mat';
%SessionName = 'S3/S3_20230327_r0_processed.mat'; %%session is problematic
%SessionName = 'S6/S6_20230710_r0_processed.mat';
%SessionName = 'S7/S7_20230608_r0_processed.mat';
%SessionName = 'S11/S11_20230801_r0_processed.mat';
%SessionName = 'S12/S12_20230731_r0_processed.mat';

%APC
%SessionName = 'Q3/Q3_20221019_r0_processed.mat'; %file corrupted
%SessionName = 'Q4/Q4_20221109_r0_processed.mat';
%SessionName = 'Q8/Q8_20221207_r0_processed.mat';
%SessionName = 'Q9/Q9_20221116_r0_processed.mat';


%% Path
if strcmp(computer,  'MACI64')
    datapath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionName);

%% Might need to change later
handles.PlotSelectTrials.Value = 1;
handles.SortReplay.Value = 1;

%% get the processed data loaded

[TracesOut, ColNames, handles.TrialInfo, handles.SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop, handles.Tuning] = ...
    LoadProcessedDataSession(MySession); 

Nb_unit = size(handles.SingleUnits,2);
ChosenUnits = 1:Nb_unit;

%% get the closed loop tuning curve
%TuningCurve = location x (mean std) x units x odor > for PG code
[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves_MD(SessionName, ChosenUnits, 'tuningbins', 15);
% in the MD version of the GetOdorTuningCurves also contains the median and sample number (idx) 
% which are useful for CI95 calculation
% dim are: location x (mean median std idx) x units x odor 

%% check that its a halt session - if not disable Halt only options
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')) || ...
            any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        handles.OnlyHaltRelated.Value = 1; 
        handles.OnlyHaltRelated.Enable = 'on';
        handles.OdorList = mode(...
            [ handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip')); ...
            handles.TrialInfo.Odor(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))] ...
            );
else
    handles.OnlyHaltRelated.Value = 0; 
    handles.OnlyHaltRelated.Enable = 'off'; 
    handles.OdorList = [1 2 3];
end

%% Get Halt Odor and Location
whichodor = handles.OdorList;

perturbationTrials = intersect(find(strncmpi(handles.TrialInfo.Perturbation(:,1),'Halt-Flip',9)), ...
find(handles.TrialInfo.Odor==whichodor));
haltlocation = handles.TrialInfo.Perturbation{perturbationTrials(2),2}(3);
%% Closed Loop: Get all spikes, all units aligned to trials  

[handles.AlignedSniffs, handles.sniffAlignedSpikes, handles.trialAlignedSpikes, ...
    handles.whichtetrode, handles.Events, handles.EventsPhase, handles.TrialInfo] = ...
    TrialAndSniffAlignedSpikeTimes(handles.SingleUnits,TTLs,size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession);

%% same for replays
if ~isempty(OpenLoop)
    
    % sniffs
    [handles.ReplayAlignedSniffs, handles.SniffAlignedReplaySpikes, handles.ReplayInfo] = ...
        SniffAlignedSpikeTimes_Replays(handles.SingleUnits,TTLs,ReplayTTLs,handles.TrialInfo,OpenLoop,MySession);
    
    % regular trials
    if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
        [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayTrialInfo] = ...
            PerturbationReplayAlignedSpikeTimes_v2(handles.SingleUnits,TTLs,...
            ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop,MySession,'sniffwarpmethod',0);
    end
else
    handles.ReplayAlignedSniffs = [];
end

%% also for passive tuning
handles.TuningSniffs = PassiveTuningSniffs(handles.Tuning,MySession);

%% pseudorandomtuning trials
if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
    [handles.PseudoRandomTuningSpikes] = ...
        TrialAlignedSpikeTimes_Tuning(handles.SingleUnits,handles.Tuning.TTLs);
    
    % transition markers
    odorTS(1,1) = handles.Tuning.extras.sessionsettings(1,4); % w.r.t. trial start (pre-odor)
    nLocations = size(handles.Tuning.extras.sequence,2) - 2;
    LocationShifts = 0; 
    for i = 1:nLocations
        if i == 1
            LocationShifts(i,1) = -handles.Tuning.extras.sessionsettings(1,3); % w.r.t. trial start (settle)
            LocationShifts(i,2) = sum(handles.Tuning.extras.sessionsettings(1,[4,5])); % w.r.t. trial start (pre + odor)
        else
            LocationShifts(i,1) = LocationShifts(i-1,2);
            LocationShifts(i,2) = LocationShifts(i,1) + sum(handles.Tuning.extras.sessionsettings(1,[3,5])); % settle + odor
        end
    end
    odorTS(1,2) = LocationShifts(end,2);
    handles.TuningTiming.LocationShifts = LocationShifts/1000; % in s
    handles.TuningTiming.Odor = odorTS/1000; % in s

end

%% Get FR aligned to perturbation start for all types of trials (CL, Halt, CL replay, Halt replay
%And FR aligned to halt location for tuning

for whichUnit = 1:Nb_unit
    for i = 1:numel(handles.OdorList)
        
        %define some variables
        whichodor = handles.OdorList(i);
        AlignType = 6; %perturbation start
        mywin = 1000; % in ms
        stepsize = 1000/SampleRate;


        %% Baseline trials (CL and halt)
        [trialsdone, AlignedFRs, BinOffset, AlignedPerturbationFRs, RawSpikeCounts] = PlotFullSession(-whichUnit, whichodor, handles.trialAlignedSpikes, handles.Events, ...
            handles.TrialInfo, handles.TrialInfo.InZone, AlignType, 'plotspikes', 0, ...
            'trialfilter', handles.PlotSelectTrials.Value, 'psth',1);

        % count spikes in a 1000ms window after perturbation start
        bins = [stepsize:stepsize:mywin] - BinOffset;
        bins = bins/stepsize;
        foo = round(numel(bins)/2);
        if size(AlignedPerturbationFRs,2)<bins(end)
            AlignedPerturbationFRs = horzcat(AlignedPerturbationFRs,...
                zeros(size(AlignedPerturbationFRs,1),(bins(end)- size(AlignedPerturbationFRs,2))));
        end
        for tz = 1:12
            % get odor response during 'windowsize' after odor ON
            AreaUnderCurve.Odor(tz,whichUnit) = sum(AlignedFRs(tz,bins)); %TZ x unit
            AreaUnderCurve.Halt(tz,whichUnit) = sum(AlignedPerturbationFRs(tz,bins)); %TZ x unit

        end

        AreaUnderCurve.Odor(:,whichUnit) = AreaUnderCurve.Odor(:,whichUnit)/(mywin/stepsize);
        AreaUnderCurve.Halt(:,whichUnit) = AreaUnderCurve.Halt(:,whichUnit)/(mywin/stepsize);

        %% passive CL replay and passive halts 
        if ~isempty(handles.ReplayAlignedSniffs)
            [perturbationreplaysadded, PassiveReplayFRs, PerturbationReplayFRs, BinOffset] = AddPerturbationReplay2FullSession_v2(trialsdone, -whichUnit, whichodor, handles.ReplayAlignedSpikes, ...
                handles.ReplayEvents, handles.ReplayTrialInfo, handles.ReplayTrialInfo.InZone, AlignType, handles.SortReplay.Value, ...
                'trialfilter', handles.PlotSelectTrials.Value, 'plotspikes', 0, 'psth',1);

            trialsdone = trialsdone + perturbationreplaysadded;
        

        %count spikes in a 1000ms window after perturbation start
        bins = [stepsize:stepsize:mywin] - BinOffset;
        bins = bins/stepsize;
        foo = round(numel(bins)/2);
        %hack since PerturbationReplayFRs doesn't go up to tz 12
        if size(PerturbationReplayFRs,2) <12
            PerturbationReplayFRs{12} = [];
        end
        if size(PassiveReplayFRs,2) <12
            PassiveReplayFRs{12} = [];
        end
        %if FRs are too short
        for tz = 1:12
            if isempty(PerturbationReplayFRs{tz})
                PerturbationReplayFRs{tz} = zeros(bins(end),1);
            end
            if isempty(PassiveReplayFRs{tz})  
                PassiveReplayFRs{tz}= zeros(bins(end),1); 
            end
            if size(PerturbationReplayFRs{tz},1)<bins(end)
                PerturbationReplayFRs{tz} = vertcat(PerturbationReplayFRs{tz},...
                     zeros((bins(end)- size(PerturbationReplayFRs{tz},1)),1));
            end
            if size(PassiveReplayFRs{tz},1)<bins(end)
                PassiveReplayFRs{tz} = vertcat(PassiveReplayFRs{tz},...
                     zeros((bins(end)- size(PassiveReplayFRs{tz},1)),1));
                   
            end

        end

        for tz = 1:12
            % get odor response during 'windowsize' after odor ON
            AreaUnderCurve.PassiveOdor(tz,whichUnit) = sum(PassiveReplayFRs{tz}(bins)); %TZ x unit
            AreaUnderCurve.PassiveHalt(tz,whichUnit) = sum(PerturbationReplayFRs{tz}(bins)); %TZ x unit
        end

        AreaUnderCurve.PassiveOdor(:,whichUnit) = AreaUnderCurve.PassiveOdor(:,whichUnit)/(mywin/stepsize);
        AreaUnderCurve.PassiveHalt(:,whichUnit) = AreaUnderCurve.PassiveHalt(:,whichUnit)/(mywin/stepsize);

        end
        %% add tuning trials 
        if any(handles.Tuning.extras.sequence(:,1)==800) % pseudorandom tuning
            TuningAlignType = 1000 + haltlocation;
            LocationDuration = mode(diff(handles.TuningTiming.LocationShifts'));
            [trialsdone, FR, BinOffset] = PlotRandomTuningTrials(trialsdone, -whichUnit, whichodor, handles.PseudoRandomTuningSpikes, ...
                handles.TuningTiming, handles.Tuning.extras.sequence, TuningAlignType, LocationDuration, [0 1], 'plotspikes', 0, 'psth',1);
        else
            [trialsdone, FR, BinOffset] = PlotTuningTrials(trialsdone, -whichUnit, whichodor, handles.SingleUnits, handles.Tuning.TTLs, ...
                'plotspikes', 0, 'selectlocation', haltlocation,'psth',1);
        
        end
        %count spikes in a 1000ms window
        bins = [stepsize:stepsize:mywin] - BinOffset;
        bins = bins/stepsize;
        foo = round(numel(bins)/2);
         if size(FR,2)<bins(end) 
             FR = horzcat(FR,zeros(size(FR,1),(bins(end)- size(FR,2))));
         end
         if isempty(FR)
             FR = zeros(1,bins(end));
         end
         % get odor response during 'windowsize' after odor ON
        AreaUnderCurve.Tuning(whichUnit) = sum(FR(bins)); %unit
        AreaUnderCurve.Tuning(:,whichUnit) = AreaUnderCurve.Tuning(:,whichUnit)/(mywin/stepsize);
    end
end



%% COMPARE MEAN PERTURBATION RESPONSE TO MIRROR ODOR LOCATION

for i = 1:Nb_unit
    % consider comparing confidance intervals of the median instead of
    % the mean for non normal distributions
    % mean includes outliers more thna median 
    % but we also have very little samples for Odor response 
    nsamps_halt = length(AreaUnderCurve.Halt([1 5 9],i));
    mean_halt = mean(AreaUnderCurve.Halt([1 5 9],i));
    median_halt = median(AreaUnderCurve.Halt([1 5 9],i));
    SEM_halt = std(AreaUnderCurve.Halt([1 5 9],i))/sqrt(nsamps_halt);
    ts_halt = tinv([0.025  0.975],(nsamps_halt-1));
    CI95_halt = mean_halt + ts_halt*SEM_halt; 
    %CI95_halt = median_halt + ts_halt*SEM_halt; 
    
    nsamps_odor = length(AreaUnderCurve.Odor(11:12,i));
    mean_odor = mean(AreaUnderCurve.Odor(11:12,i));
    median_odor = median(AreaUnderCurve.Odor(11:12,i));
    SEM_odor = std(AreaUnderCurve.Odor(11:12,i))/sqrt(nsamps_odor);
    ts_odor = tinv([0.025  0.975],(nsamps_odor-1));
    CI95_odor = mean_odor + ts_odor*SEM_odor; 
    %CI95_odor = median_odor + ts_odor*SEM_odor; 

    if  CI95_halt(2) < CI95_odor(1)
        mirror_modulation(i) = 2;

    elseif CI95_halt(1) > CI95_odor(2)
        mirror_modulation(i) = 1;
    else
        mirror_modulation(i) = 0;
    end
    
end

%% COMPARE MEAN PERTURBATION RESPONSE TO TUNING CURVE LOCATION
haltlocation = 30;
WhichBin = find(mean(XBins,2)==haltlocation);

for i = 1:Nb_unit
    nsamps_halt = length(AreaUnderCurve.Halt([1 5 9],i));
    mean_halt = mean(AreaUnderCurve.Halt([1 5 9],i));
    median_halt = median(AreaUnderCurve.Halt([1 5 9],i));
    SEM_halt = std(AreaUnderCurve.Halt([1 5 9],i))/sqrt(nsamps_halt);
    ts_halt = tinv([0.025  0.975],(nsamps_halt-1));
    CI95_halt = mean_halt + ts_halt*SEM_halt; 
    %CI95_halt = median_halt + ts_halt*SEM_halt; 

    nsamps_tuning = squeeze(TuningCurve.ClosedLoopFull(WhichBin,4,i,whichodor));
    mean_tuning = squeeze(TuningCurve.ClosedLoopFull(WhichBin,1,i,whichodor));
    median_tuning = squeeze(TuningCurve.ClosedLoopFull(WhichBin,2,i,whichodor));
    std_tuning = squeeze(TuningCurve.ClosedLoopFull(WhichBin,3,i,whichodor));
    SEM_tuning = std_tuning/sqrt(nsamps_tuning);
    ts_tuning = tinv([0.025  0.975],(nsamps_tuning-1));
    CI95_tuning = mean_tuning + ts_tuning*SEM_tuning; 
    %CI95_tuning = median_tuning + ts_tuning*SEM_tuning; 
    
    if CI95_halt(2) < CI95_tuning(1)
        location_modulation(i) = 2;
    elseif CI95_halt(1) > CI95_tuning(2)
        location_modulation(i) = 1;
    else
        location_modulation(i) = 0;
    end
end

%% COMPARE MEAN PERTURBATION RESPONSE TO PERTURBATION REPLAY - need to change to have all repeats
if ~isempty(handles.ReplayAlignedSniffs)
for i = 1:Nb_unit
    nsamps_halt = length(AreaUnderCurve.Halt([1 5 9],i));
    mean_halt = mean(AreaUnderCurve.Halt([1 5 9],i));
    median_halt = median(AreaUnderCurve.Halt([1 5 9],i));
    SEM_halt = std(AreaUnderCurve.Halt([1 5 9],i))/sqrt(nsamps_halt);
    ts_halt = tinv([0.025  0.975],(nsamps_halt-1));
    CI95_halt = mean_halt + ts_halt*SEM_halt; 
    %CI95_halt = median_halt + ts_halt*SEM_halt; 

    nsamps_halt_p = length(AreaUnderCurve.PassiveHalt([1 5 9],i));
    mean_halt_p = mean(AreaUnderCurve.PassiveHalt([1 5 9],i));
    median_halt_p = median(AreaUnderCurve.PassiveHalt([1 5 9],i));
    SEM_halt_p = std(AreaUnderCurve.PassiveHalt([1 5 9],i))/sqrt(nsamps_halt_p);
    ts_halt_p = tinv([0.025  0.975],(nsamps_halt_p-1));
    CI95_halt_p = mean_halt_p + ts_halt_p*SEM_halt_p; 
    %CI95_halt = median_halt + ts_halt*SEM_halt; 
    


    
    if CI95_halt(2) < CI95_halt_p(1)
        replay_modulation(i) = 2;
    elseif CI95_halt(1) > CI95_halt_p(2)
        replay_modulation(i) = 1;
    else
        replay_modulation(i) = 0;
    end
end
end
%% COMPARE MEAN PERTURBATION RESPONSE TO PASSIVE TUNING LOCATION - need to change to have all repeats
haltlocation = 30;
WhichBin = find(mean(XBins,2)==haltlocation);

for i = 1:Nb_unit
    nsamps_halt = length(AreaUnderCurve.Halt([1 5 9],i));
    mean_halt = mean(AreaUnderCurve.Halt([1 5 9],i));
    median_halt = median(AreaUnderCurve.Halt([1 5 9],i));
    SEM_halt = std(AreaUnderCurve.Halt([1 5 9],i))/sqrt(nsamps_halt);
    ts_halt = tinv([0.025  0.975],(nsamps_halt-1));
    CI95_halt = mean_halt + ts_halt*SEM_halt; 
    %CI95_halt = median_halt + ts_halt*SEM_halt; 

    
    mean_tuning = AreaUnderCurve.Tuning(i);

    
    if CI95_halt(2) < mean_tuning
        tuning_modulation(i) = 2;
    elseif CI95_halt(1) > mean_tuning
        tuning_modulation(i) = 1;
    else
        tuning_modulation(i) = 0;
    end
end

%% % PLOTTING %%%

%%  PLOTTING FOR MIRROR AND TUNING CURVE COMPARISONS

figure(1) % mean perturbation response vs. mean mirror location response
title('halt vs mirror location response')
for i = 1:Nb_unit
    plot(mean(AreaUnderCurve.Odor(11:12,i)),mean(AreaUnderCurve.Halt([1 5 9],i)),...
        'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
    set(gca,'XLim',[0 40],'YLim',[0 40])
    line([0 40],[0 40],'Color','k')
    hold on
end


figure(2) % expected (based on location tuning curve vs actual 
title('halt vs tuning curve response')
haltlocation = 30;
WhichBin = find(mean(XBins,2)==haltlocation);
Expected = squeeze(TuningCurve.ClosedLoopFull(WhichBin,:,:,whichodor))';
Actual = squeeze(mean(AreaUnderCurve.Halt([1 5 9],:),1));

hold on; axis square
plot(Expected(:,1),Actual, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
set(gca,'XLim',[0 40],'YLim',[0 40])
line([0 40],[0 40],'Color','k')


%% Comparison to passive halts
if ~isempty(handles.ReplayAlignedSniffs)
figure(3) %mean active perturbation response vs mean passive perturbation response
title('active halt vs passive halt response')
hold on; axis square
plot(mean(AreaUnderCurve.PassiveHalt([1 5 9],:)),mean(AreaUnderCurve.Halt([1 5 9],:)),...
    'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
set(gca,'XLim',[0 40],'YLim',[0 40])
line([0 40],[0 40],'Color','k')
end 

figure(4) %mean perturbation response vs mean pseudorandom tuning response at same location
title('halt vs pseudorandom tuning response')
hold on; axis square
plot(AreaUnderCurve.Tuning(:),mean(AreaUnderCurve.Halt([1 5 9],:),1),...
    'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
set(gca,'XLim',[0 40],'YLim',[0 40])
line([0 40],[0 40],'Color','k')


%% Comparison between location tuning curve to pseudorandom tuning
figure(5) %mean perturbation response vs mean tuning response at same location
title('location tuning curve vs pseudorandom tuning response')

plot(Expected,mean(AreaUnderCurve.Halt([1 5 9],:),1),...
    'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
set(gca,'XLim',[0 40],'YLim',[0 40])
line([0 40],[0 40],'Color','k')
