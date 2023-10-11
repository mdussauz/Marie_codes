function [trialsdone, AlignedFRs, AlignedPerturbationFRs, RawSpikeCounts, RawPerturbationSpikeCounts] = ...
    TrialAlignedActivity(whichUnit, whichodor, AlignedSpikes, Events, TrialInfo, AlignTo, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;

params.addParameter('sniffaligned', false, @(x) islogical(x) || x==0 || x==1);
params.addParameter('sniffscalar', 3, @(x) isnumeric(x));


% extract values from the inputParser
params.parse(varargin{:});
sniffaligned = params.Results.sniffaligned;
sniffscalar = params.Results.sniffscalar;

whichUnit = abs(whichUnit);

%% hack to prevent OL-Template trials to be considered as perturbed trials
f = find(strcmp(TrialInfo.Perturbation(:,1),'OL-Template'));
if ~isempty(f)
    for i = 1:numel(f)
        TrialInfo.Perturbation{f(i),1} = [];
    end
end

%%
thisUnitSpikes = AlignedSpikes(:,whichUnit);

%% get the trial sorting order
whichTrials = intersect(find(cellfun(@isempty, TrialInfo.Perturbation(:,1))), ...
    find(TrialInfo.Odor==whichodor));
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) TrialInfo.Duration(whichTrials)]; 
whichTrials = sortrows(whichTrials,2);

for tz = 1:12
    whichTrials(whichTrials(:,2)==tz,:) = sortrows(whichTrials(whichTrials(:,2)==tz,:),3);
end

% also collect perturbation trials
perturbationTrials = intersect(find(~cellfun(@isempty, TrialInfo.Perturbation)), ...
    find(TrialInfo.Odor==whichodor));
perturbationTrials = intersect(find(~strcmp(TrialInfo.Perturbation(:,1),'OL-Replay')), ...
    perturbationTrials);
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

allTrials = vertcat(whichTrials, perturbationTrials);

%% Events to align to
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Xlims = [-1.2 -1];
        Offset = 0*myEvents(:,1);
    case 2 % odor ON
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
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
        % offset all events with ON timestamp
        myEvents = myEvents - Offset;
        Xlims = [-1.2 -1] - 1;
    case 6 % perturbation start
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

%%
x = size(whichTrials,1);
% calculate PSTH
AlignedFRs = []; RawSpikeCounts = [];
BinOffset = round(Xlims(1)*1000);

for TZ = 1:12
    thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
    Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
    [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
    AlignedFRs(TZ,1:numel(myFR)) = myFR;
    RawSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
end

entries_done = TZ;

%%
% calculate PSTH
AlignedPerturbationFRs = [];
for TZ = 1:12
    thisTZspikes = thisUnitSpikes(perturbationTrials(find(perturbationTrials(:,2)==TZ),1));
    Events2Align = Offset(x+find(perturbationTrials(:,2)==TZ),1);
    [myFR, myPSTH] = MakePSTH_v3(thisTZspikes,Events2Align,BinOffset,'downsample',500);
    AlignedPerturbationFRs(TZ,1:numel(myFR)) = myFR;
    RawPerturbationSpikeCounts(TZ,1:numel(myPSTH)) = myPSTH;
end

trialsdone = size(allTrials,1);
end 

