function [allclusters] = MakeSpikeStucture(SessionInfo,BrainRegion, Mouse)
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

if strcmp(computer, 'PCWIN64')
    %addpath(genpath('Z:\mdussauz\PhotonCerber_Stimuli_on_server'))
    addpath(genpath('Z:\mdussauz\PhotonCerber_Stimuli_on_server'))
    addpath(genpath('C:\Users\Marie\Documents\data\PhotonCerber_Stimuli_on_server'))
    
    KSpath = 'Z:\mdussauz\ephysdata\Conc_id_exp'; 
    KSsortedpath = 'Z:\mdussauz\ConcId\Sorted'; 
else
    addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
end
%% Get Trial Timestamps from Open Ephys Events file

allclusters = struct('id',[],'spikecount',[],'spikes',[], 'stimulus', []); %initiate
for i = 1: length(myKsDir)
    i
    stimfilename{i}
    [StimTime, StimList, Nrepeats] = ReadStimFile([stimfilename{i}, '.txt']);
    NOdors = numel(unique(StimList));
    NTrials= numel(StimList);
    
    filenameKS = fullfile(KSpath, myKsDir{i},'Record Node 106');
    filenameKS
    %[data, timestamps, info] = load_open_ephys_data(filenameKS); % data has channel IDs
    [data, timestamps, info] = load_open_ephys_data(fullfile(filenameKS,'all_channels.events')); % data has channel IDs
    
    % adjust for clock offset between open ephys and kilosort
    %[offset] = AdjustClockOffset(myKsDir{i});
    [offset] = AdjustClockOffset(filenameKS);
    %offset = offset/sampleRate;
    
    
    %% Get Events and correct for ephys offset
    [Events] = ParseOpenEphysEvents(fullfile(filenameKS,'all_channels.events'));
    % trial events are in Channel0, odor events are in Channel1

    TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
    % delete first entry - just an empty trigger
    TrialTimestamps(1:2,:) = [];
    
    OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
    OdorTimestamps(1,:) = [];
    
    % should i correct for offset???
      TrialTimestamps = TrialTimestamps - offset; %adjust for offset
      OdorTimestamps = OdorTimestamps - offset; % adjust for offset
    
%     %% In new experiment there seems to be an extra trial and odor on/off
%     if OdorTimestamps(1,2) - OdorTimestamps(1,1) < 2
%         TrialTimestamps(1,:) = [];
%         OdorTimestamps(1,:) = [];
%     end
        
    %% Load data from kilosort/phy
    filenameKSsorted = fullfile(KSsortedpath, myKsDir{i}); 
    %sp = loadKSdir(myKsDir{i});
    sp = loadKSdir(filenameKSsorted);
    
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
                tstart = TrialTimestamps(mytrial,1)-15;
                tstop = TrialTimestamps(mytrial,2)+15;
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
