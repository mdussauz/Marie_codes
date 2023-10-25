function [trialsdone, AlignedFRs, AlignedPerturbationFRs, RawSpikeCounts, RawPerturbationSpikeCounts] = ...
    EventAlignedActivity(whichUnit, whichodor, AlignedSpikes, Events, TrialInfo, AlignTo, window, varargin)

narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;

params.addParameter('sniffscalar', 3, @(x) isnumeric(x));
params.addParameter('TZsorted', 0, @(x) isnumeric(x)); 


% extract values from the inputParser
params.parse(varargin{:});
sniffscalar = params.Results.sniffscalar;
TZsorted = params.Results.TZsorted;

if sniffscalar~=0
    sniffaligned = 1;
else
    sniffaligned = 0;
end 

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
perturbationTrials = [perturbationTrials TrialInfo.TargetZoneType(perturbationTrials) TrialInfo.Duration(perturbationTrials)]; 
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

%% Get offset to align spikes to each event
myEvents = Events(allTrials(:,1),:);
switch AlignTo
    case 1 % to trial ON
        Offset = 0*myEvents(:,1);
    case 2 % odor ON
        if ~sniffaligned
            odorON = myEvents(:,1);
        else
            odorON = floor(myEvents(:,1));
        end
        Offset = odorON;
    case 3 % trial OFF
        if ~sniffaligned
            TrialOFF = myEvents(:,3);
        else
            TrialOFF = floor(myEvents(:,3));
        end
        Offset = TrialOFF;
    case 4 % reward
        if ~sniffaligned
            Reward = myEvents(:,3);
        else
            Reward = floor(myEvents(:,3));
        end
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
end 

% perturbation start - for now not being used
if ~isempty(perturbationTrials) && AlignTo == 6  %if there are perturbation trials
    if ~sniffaligned
        Offset = myEvents(:,5);
    else
        Offset = floor(myEvents(:,5));
    end
    % offset all events with ON timestamp
end

if sniffaligned
    window = sniffscalar*window;
end


if TZsorted 


%% calculate PSTH for each TZ for normal trials 
% initialize
BinOffset = window(1); 
windowsize = length(window(1):window(2)); 
AlignedFRs = zeros(12,windowsize); RawSpikeCounts = zeros(12,windowsize);
myPSTH = zeros(1,windowsize);  myFR = zeros(1,windowsize);

for TZ = 1:12
    thisTZspikes = thisUnitSpikes(whichTrials(find(whichTrials(:,2)==TZ),1));
    if ~isempty(thisTZspikes) % sometimes no spikes in trial
        Events2Align = Offset(find(whichTrials(:,2)==TZ),1);
        [myFR, myPSTH] = MakePSTH_MD(thisTZspikes,Events2Align,BinOffset,window,'kernelsize',25);

        AlignedFRs(TZ,:) = myFR;
        RawSpikeCounts(TZ,:) = myPSTH;
    end

end


%% calculate PSTH for each TZ for perturbed trials
x = size(whichTrials,1); %nb of normal trials
AlignedPerturbationFRs = zeros(12,windowsize); RawPerturbationSpikeCounts = zeros(12,windowsize);
myPSTHp = zeros(1,windowsize);  myFRp = zeros(1,windowsize);

if ~isempty(perturbationTrials) %if there are perturbation trials

    for TZ = 1:12
        thisTZspikes = thisUnitSpikes(perturbationTrials(find(perturbationTrials(:,2)==TZ),1));
        if ~isempty(thisTZspikes) % sometimes no spikes in trial
            Events2Align = Offset(x+find(perturbationTrials(:,2)==TZ),1);
            [myFRp, myPSTHp] = MakePSTH_MD(thisTZspikes,Events2Align,BinOffset,window,'kernelsize',25);

            AlignedPerturbationFRs(TZ,:) = myFRp;
            RawPerturbationSpikeCounts(TZ,:) = myPSTHp;
        end
    end


end

else 

    %% calculate PSTH for each trial normal trials 
% initialize
BinOffset = window(1); 
windowsize = length(window(1):window(2)); 
x = size(whichTrials,1); %nb of normal trials
AlignedFRs = zeros(x,windowsize); RawSpikeCounts = zeros(x,windowsize);
myPSTH = zeros(1,windowsize);  myFR = zeros(1,windowsize);

for thisTrial = 1:x
    thisTrialSpikes = thisUnitSpikes{whichTrials(thisTrial,1)}{1}; 
    if ~isempty(thisTrialSpikes) % sometimes no spikes in trial
    Events2Align = Offset(thisTrial);
    [myFR, myPSTH] = MakePSTH_MD(thisTrialSpikes,Events2Align,BinOffset,window,'kernelsize',25);

    AlignedFRs(thisTrial,:) = myFR;
    RawSpikeCounts(thisTrial,:) = myPSTH;
    end

end


%% calculate PSTH for each perturbed trials
x = size(whichTrials,1); %nb of normal trials
y = size(perturbationTrials,1); %nb of perturbed trials
AlignedPerturbationFRs = zeros(y,windowsize); RawPerturbationSpikeCounts = zeros(y,windowsize);
myPSTHp = zeros(1,windowsize);  myFRp = zeros(1,windowsize);

if ~isempty(perturbationTrials) %if there are perturbation trials

    for thisPertubTrial = 1:y
        thisTrialSpikes = thisUnitSpikes{perturbationTrials(thisPertubTrial,1)}{1}; 
        if ~isempty(thisTrialSpikes) % sometimes no spikes in trial
        Events2Align = Offset(x+thisPerturbTrial);
        [myFRp, myPSTHp] = MakePSTH_MD(thisTrialspikes,Events2Align,BinOffset,window,'kernelsize',25);

        AlignedPerturbationFRs(thisPertubTrial,:) = myFRp;
        RawPerturbationSpikeCounts(thisPertubTrial,:) = myPSTHp;
        end
    end

    
end

end 

trialsdone = size(allTrials,1);
end 

