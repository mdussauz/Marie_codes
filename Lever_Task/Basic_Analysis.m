%Basic_Analysis 

%% Select a processed file to analyze

% [WhichSession, SessionPath] = uigetfile(); % to be used later for more flexibility

WhichSession = 'O3_20211005_r0_processed';
SessionPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior\O3';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);

%% get the data loaded
MySession = handles.WhereSession.String;
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession);  % loads relevant variables

% SingleUnits contains the spiketimes of all units from start of trial up
% to the start of next trial 

% startoffset = Time window (in seconds) preceding trial start used for extracting Traces
startoffset = 1;

%% Get all spikes, all units aligned to trials - both for closed-loop and replays
[handles.AlignedSpikes, handles.Events, handles.whichtetrode] = TrialAlignedSpikeTimes(SingleUnits,TTLs,...
    size(handles.TrialInfo.TrialID,2),handles.TrialInfo,MySession); 

% trialalignedspikes = spiketimes relative to the most recent trial''s'' start timestamp');

% Perturbations - Might be useful at some point but not used right now
if any(strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        ReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events);
end

if any(strcmp(handles.TrialInfo.Perturbation(:,1),'Halt-Flip-Template'))
    [handles.ReplayAlignedSpikes, handles.ReplayEvents, handles.ReplayInfo] = ...
        PerturbationReplayAlignedSpikeTimes(SingleUnits,TTLs,...
        ReplayTTLs,handles.TrialInfo,handles.Events,OpenLoop);
end

%% Getting number of units in session > will be useful to normalize results .i.e responsiveness 
handles.NumUnits.String = num2str(size(SingleUnits,2));

N = handles.NumUnits.String; % renaming var
MyUnits = (1:N);

%% Getting the FR for different windows

Binsize  = 50; 

for whichodor = 1:3 % for every odor 
    % 1. Pick the right trials - no perturbation, and a given odor
    whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation)), ...
        find(TrialInfo.Odor==whichodor));
    XVar = []; YVar = [];
    for whichUnit = 1:numel(MyUnits) % for every Unit
        counts = [1 0];
        thisUnitSpikes = AlignedSpikes(:,MyUnits(whichUnit));
        for x = 1:size(whichTrials,1) % every trial
            % only keep the PSTH and behavioral variables from odor start until Trial OFF
            t1 = round((startoffset + Events(whichTrials(x,1),1))*SampleRate);
            t2 = round(Events(whichTrials(x,1),3)*SampleRate);
            if mod(numel(t1:t2),Binsize/2)
                t2 = t1 + (Binsize/2)*floor(numel(t1:t2)/(Binsize/2)) - 1;
            end
            newSamples = numel(t1:t2)/(Binsize/2);
            counts(2) = counts(2) + newSamples; % one extra bin for pre-odor start spikes
            if whichUnit == 1
                % get behavioral variable
                myMotor = mean(reshape(Traces.Motor{whichTrials(x,1)}(t1:t2),Binsize/2,[]));
                XVar(counts(1):counts(2),1) = myMotor';
            end
            thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1}; % aligned to trial ON
            thisTrialSpikes = thisTrialSpikes - Events(whichTrials(x,1)); % subtract odor start - odor start becomes zero
            % convert to ms
            thisTrialSpikes = ceil(thisTrialSpikes*1000/Binsize); % at binned resolution
            thisTrialSpikes(thisTrialSpikes<1) = [];
            [C,~,ic] = unique(thisTrialSpikes);
            bin_counts = accumarray(ic,1);
            myRaster = [];
            if ~isempty(C)
                myRaster(C) = bin_counts;
            else
                myRaster = zeros(1,newSamples);
            end
            if numel(myRaster)<diff(counts)+1
                myRaster = [myRaster, zeros(1,(diff(counts) + 1 - numel(myRaster)))];
            end
            YVar(counts(1):counts(2),whichUnit) = myRaster(1:(1+diff(counts)));
            
            counts(1) = counts(2)+1;
        end
    end
    [Curve_CL{whichodor}, XBins] = SmellocatorTuning('Odor',125-XVar,1000*(YVar/Binsize));
end



%% Information I need here 

whichUnit = 6; %  SPECIFY 
thisUnitSpikes = handles.AlignedSpikes(:,whichUnit);

whichodor = 3; % SPECIFY

% get the trial sorting order
whichTrials = intersect(find(cellfun(@isempty, handles.TrialInfo.Perturbation(:,1))), ...
    find(handles.TrialInfo.Odor==whichodor));
whichTrials = [whichTrials handles.TrialInfo.TargetZoneType(whichTrials) handles.TrialInfo.Duration(whichTrials)];
whichTrials = sortrows(whichTrials,2);

for tz = 1:12
    whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
end

% also collect perturbation trials
perturbationTrials = intersect(find(~cellfun(@isempty, handles.TrialInfo.Perturbation)), ...
    find(handles.TrialInfo.Odor==whichodor));
perturbationTrials = intersect(find(~strcmp(handles.TrialInfo.Perturbation(:,1),'OL-Replay')), ...
    perturbationTrials);
perturbationTrials = [perturbationTrials handles.TrialInfo.TargetZoneType(perturbationTrials) handles.TrialInfo.Duration(perturbationTrials)]; %#ok<AGROW>
perturbationTrials = sortrows(perturbationTrials,2);
for tz = 1:12
    perturbationTrials(perturbationTrials(:,2)==tz,:) = sortrows(perturbationTrials(perturbationTrials(:,2)==tz,:),3);
end

allTrials = vertcat(whichTrials, perturbationTrials);

for AlignTo = 1:5 % for now ignoring perturbations % need to add a clause to check if there are perturbations
    % Plot all events
    myEvents = handles.Events(allTrials(:,1),:);
    switch AlignTo
        case 1 % to trial ON
            Xlims = [-1.2 -1];
            Offset = 0*myEvents(:,1);
        case 2 % odor ON
            odorON = myEvents(:,1);
            myEvents(:,1) = 0; % replace odorON with TrialON
            % offset all events with ON timestamp
            myEvents = myEvents - odorON;
            Xlims = [-1.2 -1];
            Offset = odorON;
        case 3 % trial OFF
            TrialOFF = myEvents(:,3);
            myEvents(:,3) = 0; % replace TrialOFF with TrialON
            % offset all events with ON timestamp
            myEvents = myEvents - TrialOFF;
            Xlims = [-1.2 -1] - 4;
            Offset = TrialOFF;
        case 4 % reward
            Reward = myEvents(:,3);
            myEvents(:,2) = 0; % replace Reward with TrialON
            % offset all events with ON timestamp
            myEvents = myEvents - Reward;
            Xlims = [-1.2 -1] - 4;
            Offset = Reward;
        case 5 % first TZ entry with stay > 100ms
            Offset = myEvents(:,4);
            % offset all events with ON timestamp
            myEvents = myEvents - Offset;
            Xlims = [-1.2 -1] - 1;
%         case 6 % perturbation start
%             Offset = myEvents(:,5);
%             % offset all events with ON timestamp
%             myEvents = myEvents - Offset;
%             Xlims = [-1.2 -1];
    end
    
    %% Calculate PSTH
    AlignedFRs = []; RawSpikeCounts = [];
    BinOffset = Xlims(1)*1000; % would be window start here
    
    for TZ = 1:12
        thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
        Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
        [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
        AlignedFRs(TZ,1:numel(myFR)) = myFR;
        RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
    end
    
    %%
%     figure(AlignTo)
%     for stim = 1:numStim
%         
%         subplot(, , )
%         
%         air = sum(odorspikes < 0)/50;
%         odor = sum(odorspikes > 0.3 & odorspikes < 4)/20;
%         bar(0, air);
%         hold on
%         bar(1,odor );
%         
%     end

end