function [allclusters] = MakeSpikeStucture(SessionInfo,BrainRegion, Mouse)
% -- written by MD
% -- create structure array containing: 
% 
% user must specify what to extract: 
% sessions coming from AON or APC mouse 
% whether to use all mice or just a specific mouse 
% whether to build it for a specific session   

AnalysisChoice = [BrainRegion,'_', Mouse];

%% defaults
%sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz - not used 

%% Stimulus File and KS directory 
% get filenames info based on user's choice  

switch AnalysisChoice
    case 'AON_all'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
        mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
    case 'APC_all'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion));
        mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion));
    case 'AON_E3'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E3',SessionInfo.MouseName));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E3',SessionInfo.MouseName));
        mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E3',SessionInfo.MouseName));
    case 'AON_E2'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E2',SessionInfo.MouseName));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E2',SessionInfo.MouseName));
        mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion)& strcmp('E2',SessionInfo.MouseName));
    case 'APC_E6'
        stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion)& strcmp('E6',SessionInfo.MouseName));
        myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion)& strcmp('E6',SessionInfo.MouseName));
        mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('APC',SessionInfo.BrainRegion)& strcmp('E6',SessionInfo.MouseName));
end

if strcmp(computer, 'PCWIN64') % Marie's remote desktop 
    addpath(genpath('Z:\mdussauz\PhotonCerber_Stimuli_on_server'))
    addpath(genpath('C:\Users\Marie\Documents\data\PhotonCerber_Stimuli_on_server'))
    
    KSpath = 'Z:\mdussauz\ephysdata\Conc_id_exp'; 
    KSsortedpath = 'Z:\mdussauz\ConcId\Sorted'; 
    
else % Marie's work linux machine
    addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))

    KSpath = '/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp';
    KSsortedpath = '/mnt/grid-hs/mdussauz/ConcId/Sorted'; 
end

%% Make spikes structure 
allclusters = struct('id',[],'spikecount',[],'spikes',[], 'stimulus', [], 'settings',[], 'repeats', []); 

for i = 1: length(myKsDir) %loop through each rec file
    [StimTime, StimList, Nrepeats] = ReadStimFile([stimfilename{i}, '.txt']);
    NTrials= numel(StimList); 
    
    % sanity check to verify whether all repeats of each odors were done
    NOdors = numel(unique(StimList));
    ExpectedNTrials = NOdors*Nrepeats;
    if ExpectedNTrials ~= NTrials
        disp(['not all repeats were done in session ', stimfilename{i}])
    end 
    
    % this next line seems unused but double check:
    %[data, timestamps, info] = load_open_ephys_data(fullfile(filenameKS,'all_channels.events')); % data has channel IDs
     
    %% Get trial and odor events and correct for ephys offset from open ephys file
    filenameKS = fullfile(KSpath,mousename{i}, myKsDir{i},'Record Node 106');
    [Events] = ParseOpenEphysEvents(fullfile(filenameKS,'all_channels.events'));

    % trial events are in Channel0: 
    TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
    
    C = strsplit(myKsDir{i},'-'); %split the KsDir name
    ExpYear = str2double(C{1}); % year of experiment
    
    if ExpYear >= 2022 % for experiments that happened on/after 2022
        % delete 1st entry - just an empty trigger - and
        % delete 2nd entry - extra trigger on/off from DAC reset (2022 change):
        TrialTimestamps(1:2,:) = [];
        
        % odor events are in Channel1:
        OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
        % delete 1st entry - extra odor on/off from DAC reset (2022 change):
        OdorTimestamps(1,:) = [];
        
    else % for experiments before 2022
        TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
    end
    
    % Sanity check - make sure number of TTLs are matching number of trials
    if ExpectedNTrials ~= length(TrialTimestamps)
        disp(['error in # of Trial TTLs in session ', stimfilename{i}])
    end 
    if ExpectedNTrials ~= length(OdorTimestamps)
        disp(['error in # of Odor TTLs in session ', stimfilename{i}])
    end
    
    % adjust for clock offset between open ephys (events) and kilosort (spike times)
    [offset] = AdjustClockOffset(filenameKS); %read 1st timestamps
    TrialTimestamps = TrialTimestamps - offset; %adjust for offset
    OdorTimestamps = OdorTimestamps - offset; % adjust for offset
            
    %% Load data from sorted session (ran through kilosort + phy)
    filenameKSsorted = fullfile(KSsortedpath, mousename{i},myKsDir{i}); 
    sp = loadKSdir(filenameKSsorted);
    
    %% Number of good units
    goodUNber = length (find(sp.cgs==2));
    disp(['found ',num2str(goodUNber),' good units']);
    
    %% Split data by clusters and by trials  
    for mycluster = 1:length(sp.cids) % loop through each cluster
        
        if sp.cgs(mycluster) == 2 % keep clusters labeled as good
            %get spike times for this cluster
            allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds
            clear goodspiketimes
            
            for mytrial = 1:NTrials % loop through each trial 
                %get trial start, stop and odor start times for this trial: 
                tstart = TrialTimestamps(mytrial,1); % optional: substract -15 to get ITI
                tstop = TrialTimestamps(mytrial,2); %optional: add +15 to get ITI
                odorstart = OdorTimestamps(mytrial,1);
                %only keep spikes in trial:
                mygoodspikes = allgoodspikes(find(allgoodspikes>=tstart & allgoodspikes<=tstop));
                mygoodspikes = mygoodspikes - odorstart; % align to every odor start
                goodspiketimes(mytrial) = {mygoodspikes}; %store spikes per trial
                %each column of goodspiketimes is a trial
                %they contain the times in sec at which spikes occured 
            end
            goodcluster(mycluster).id = sp.cids(mycluster); %store cluster id
            goodcluster(mycluster).spikecount = numel(allgoodspikes);%store number of spikes
            goodcluster(mycluster).spikes = goodspiketimes; %store spikes in trial
            goodcluster(mycluster).stimulus = StimList; %store odor info
            goodcluster(mycluster).settings = StimTime; %store trial setting info
            goodcluster(mycluster).repeats = Nrepeats; %store number of repeats
        end
    
    end
    
    allclusters =  [allclusters, goodcluster];
    allclusters(all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), allclusters, 'UniformOutput', false )), 1 )) = []; %remove empty lines
    clear goodcluster
end
