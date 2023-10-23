function [PerturbationReplayAlignedFRs,PerturbationReplayRawSpikeCounts,...
    OtherPerturbationReplayAlignedFRs,OtherPerturbationReplayRawSpikeCounts]...
    = PerturbationReplayTrialAlignedActivity(trialsdone, whichUnit,...
    whichodor, AlignedSpikes, Events, TrialInfo, AlignTo, SortTrials, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;

params.addParameter('sniffscalar', 3, @(x) isnumeric(x));

% extract values from the inputParser
params.parse(varargin{:});
sniffscalar = params.Results.sniffscalar;

if sniffscalar~=0
    sniffaligned = 1;
else
    sniffaligned = 0;
end 

%%
thisUnitSpikes = AlignedSpikes(:,whichUnit);


%% get the trial sorting order
whichTrials = intersect(find(TrialInfo.Odor==whichodor), find(~TrialInfo.Perturbed)); 
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) ...
               TrialInfo.Duration(whichTrials)]; 

% Sort trials by target zone type
if SortTrials
    for tz = 1:12
        whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
    end
end

% also collect perturbation trials
perturbationTrials = intersect(find(TrialInfo.Odor==whichodor), find(TrialInfo.Perturbed));
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) TrialInfo.Duration(perturbationTrials)]; %#ok<AGROW>
perturbationTrials = sortrows(perturbationTrials,2);
for tz = 1:12
    perturbationTrials(perturbationTrials(:,2)==tz,:) = sortrows(perturbationTrials(perturbationTrials(:,2)==tz,:),3);
end

% for offsets - sort by offset type
if any(perturbationTrials)
    if any(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template')) || ...
            any(strcmp(TrialInfo.Perturbation(:,1),'Offset-II'))
        offsetParams = cell2mat(TrialInfo.Perturbation(perturbationTrials(:,1),2));
        [~,sortidx] = sort(offsetParams(:,3));
        perturbationTrials = perturbationTrials(sortidx,:);
        offsetParams = offsetParams(sortidx,:);
        offsettypes = unique(offsetParams(:,3));
        for k = 1:numel(offsettypes)
            f = find(offsetParams(:,3)==offsettypes(k));
            perturbationTrials(f,2) = perturbationTrials(f,2) + 0.1*k;
        end
    end
end

allTrials = vertcat(perturbationTrials, whichTrials);

%% Events
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,4);
    case 2 % odor ON
        odorON = myEvents(:,1);
        if ~sniffaligned
            odorON = myEvents(:,1);
        else
            odorON = floor(myEvents(:,1));
        end
        myEvents(:,1) = 0; % replace odorON with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - odorON;
        Xlims = [-1.2 -1];
        Offset = odorON;
    case 3 % trial OFF
        TrialOFF = myEvents(:,3);
        if ~sniffaligned
            TrialOFF = myEvents(:,3);
        else
            TrialOFF = floor(myEvents(:,3));
        end
        myEvents(:,3) = 0; % replace TrialOFF with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - TrialOFF;
        Xlims = [-1.2 -1] - 4;
        Offset = TrialOFF;
    case 4 % reward
        Reward = myEvents(:,3);
        if ~sniffaligned
            Reward = myEvents(:,3);
        else
            Reward = myEvents(:,3);
        end
        myEvents(:,2) = 0; % replace Reward with TrialON
        % offset all events with ON timestamp
        myEvents = myEvents - Reward;
        Xlims = [-1.2 -1] - 4;
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        Offset = myEvents(:,4);
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
    case 6 % perturbation start
        Offset = myEvents(:,5);
        if ~sniffaligned
            Offset = myEvents(:,5);
        else
            Offset = floor(myEvents(:,5));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1];
end

if sniffaligned
    Xlims = sniffscalar*Xlims;
end

%% Calculate PSTH
% Plot Spikes
for x = 1:size(allTrials,1)

    thisTrialSpikes = thisUnitSpikes{allTrials(x,1)}{1};
    % adjust spiketimes if needed
    thisTrialSpikes = thisTrialSpikes - Offset(x);

    if x>size(perturbationTrials,1)
        % color blue on raster
        Events2Align = Offset(x);
        [myFR, myPSTH] = MakePSTH_v3(thisTrialSpikes,Events2Align,BinOffset,'downsample',500);
        PerturbationReplayAlignedFRs(x,1:numel(myFR)) = myFR;
        PerturbationReplayRawSpikeCounts(x,1:numel(myPSTH)) = myPSTH;
    else
        % color red on raster
        Events2Align = Offset(x);
        [myFR, myPSTH] = MakePSTH_v3(thisTrialSpikes,Events2Align,BinOffset,'downsample',500);
        OtherPerturbationReplayAlignedFRs(x,1:numel(myFR)) = myFR;
        OtherPerturbationReplayRawSpikeCounts(x,1:numel(myPSTH)) = myPSTH;
    end

end

end