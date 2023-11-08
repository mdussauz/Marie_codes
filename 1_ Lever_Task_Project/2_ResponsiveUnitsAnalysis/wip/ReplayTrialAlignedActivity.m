function [ActiveAlignedFRs, ActiveRawSpikeCounts, PassiveAlignedFRs, ...
    PassiveRawSpikeCounts] = ReplayTrialAlignedActivity(trialsdone, ...
    whichUnit, whichodor, AlignedSpikes, Events, TrialInfo, AlignTo, windowsize, SortTrials, varargin)

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

thisUnitSpikes = AlignedSpikes(:,whichUnit);


%% get the trial sorting order  
whichTrials = find(TrialInfo.Odor==whichodor); % both active and passive replays
whichTrials = [whichTrials TrialInfo.TargetZoneType(whichTrials) ...
               TrialInfo.Duration(whichTrials) (TrialInfo.TrialID(whichTrials)<0)']; 


%% Sort trials - first by active and passive replay
if SortTrials
whichTrials(whichTrials(:,4)==0,:) = sortrows(whichTrials(whichTrials(:,4)==0,:),2);
whichTrials(whichTrials(:,4)==1,:) = sortrows(whichTrials(whichTrials(:,4)==1,:),2);


    for tz = 1:12
        q = find((whichTrials(:,2)==tz)&(whichTrials(:,4)==0));
        whichTrials(q,:) = sortrows(whichTrials(q,:),3);
        q = find((whichTrials(:,2)==tz)&(whichTrials(:,4)==1));
        whichTrials(q,:) = sortrows(whichTrials(q,:),3);
    end
end

%% Events to align to 
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
        myEvents(:,1) = 0; % replace odorON with TrialON
        Offset = odorON;
    case 3 % trial OFF
        if ~sniffaligned
            TrialOFF = myEvents(:,3);
        else
            TrialOFF = floor(myEvents(:,3));
        end
        myEvents(:,3) = 0; % replace TrialOFF with TrialON
        Offset = TrialOFF;
    case 4 % reward
        if ~sniffaligned
            Reward = myEvents(:,3);
        else
            Reward = myEvents(:,3);
        end
        myEvents(:,2) = 0; % replace Reward with TrialON
        Offset = Reward;
    case 5 % first TZ entry with stay > 100ms
        if ~sniffaligned
            Offset = myEvents(:,4);
        else
            Offset = floor(myEvents(:,4));
        end
        % offset all events with ON timestamp
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

BinOffset = round(Xlims(1)*1000);

%% calculate PSTH
for x = 1:size(whichTrials,1)
    thisTrialSpikes = thisUnitSpikes{whichTrials(x,1)}{1};
    % adjust spiketimes if needed
    thisTrialSpikes = thisTrialSpikes - Offset(x);

    if TrialInfo.TrialID(whichTrials(x))>0
        %Color red - Active?
        Events2Align = Offset(x);
        [myFR, myPSTH] = MakePSTH_v3(thisTrialSpikes,Events2Align,BinOffset,'downsample',500);
        ActiveAlignedFRs(x,1:numel(myFR)) = myFR;
        ActiveRawSpikeCounts(x,1:numel(myPSTH)) = myPSTH;
    else
        %Color teal - Passive?
        Events2Align = Offset(x);
        [myFR, myPSTH] = MakePSTH_v3(thisTrialSpikes,Events2Align,BinOffset,'downsample',500);
        PassiveAlignedFRs(x,1:numel(myFR)) = myFR;
        PassiveRawSpikeCounts(x,1:numel(myPSTH)) = myPSTH;
    end

end


end