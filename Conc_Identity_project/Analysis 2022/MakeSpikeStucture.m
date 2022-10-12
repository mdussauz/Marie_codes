function [allcluster] = MakeSpikeStucture(BrainRegion, Mouse)
% -- written by MD
% -- create structure array for all neurons all sessions all mice
% user must specify what to extract: 
% sessions coming from AON or APC mouse 
% whether to use all mice or just a specific mouse 
% whether to build it for a specific session   

AnalysisChoice = [BrainRegion,'_', Mouse];

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Stimulus File and KS directory 

switch AnalysisChoice
    case 'AON_all'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
    case 'APC_all'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion));
    case 'AON_E3'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp( 'E3',SessionInfo.MouseName));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp( 'E3',SessionInfo.MouseName));
end

%% Get Trial Timestamps from Open Ephys Events file

allclusters = struct('id',[],'spikecount',[],'spikes',[], 'stimulus', []); %initiate
for i = 1: length(myKsDir)
    
    [StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename{i});
    NOdors = numel(unique(StimList));
    NTrials= numel(StimList);
    
    filename = fullfile(myKsDir{i},'all_channels.events');
    [data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs
    
    % adjust for clock offset between open ephys and kilosort
    [offset] = AdjustClockOffset(myKsDir{i});
    offset = offset/sampleRate;
    
    
    %% Get Events and correct for ephys offset
    [Events] = ParseOpenEphysEvents(filename);
    % trial events are in Channel0, odor events are in Channel1
    
    TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
    TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
    OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
    
    TrialTimestamps = TrialTimestamps - offset; %adjust for offset
    OdorTimestamps = OdorTimestamps - offset; % adjust for offset
    
    %% Load data from kilosort/phy
    
    sp = loadKSdir(myKsDir{i});
    
    %% Number of good units
    
    goodUNber = length (find(sp.cgs==2));
    disp(['found ',num2str(goodUNber),' good units']);
    
    %% Split data by clusters and by trials
    
    for mycluster = 1:length(sp.cids) % for each cluster
        
        if sp.cgs(mycluster) == 2 % for clusters labeled as good
            %get spike times for this cluster
            allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds
            clear goodspiketimes
            
            for mytrial = 1:NTrials
                %TrialSpikes = [];
                tstart = TrialTimestamps(mytrial,1);
                tstop = TrialTimestamps(mytrial,2);
                odorstart = OdorTimestamps(mytrial,1);
                mygoodspikes = allgoodspikes(find(allgoodspikes>=tstart & allgoodspikes<=tstop));
                mygoodspikes = mygoodspikes - odorstart; % align to every odor start!!!
                goodspiketimes(mytrial) = {mygoodspikes}; %spikes in trial
                %each column of goodspiketimes is a trial
                %for each trial gives you the spiketimes of spikes that occured
            end
            goodcluster(mycluster).id = sp.cids(mycluster);
            goodcluster(mycluster).spikecount = numel(allgoodspikes);%all spikes
            goodcluster(mycluster).spikes = goodspiketimes; %spikes in trial
            % stores goodspiketimes for each cluster before it gets cleared at
            % next iteration of the loop
            goodcluster(mycluster).stimulus = StimList;
        end
        
        
    end
    
    allclusters =  [allclusters, goodcluster];
    allclusters( all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), allclusters, 'UniformOutput', false ) ), 1 ) ) = []; %remove empty lines
    clear goodcluster
end
