function [goodcluster, muacluster] = Phy2MatProcessing(myKsDir)

%written by MD 

% processing of spike times and store them in stc based on unit types and
% trial number

% input = directory with kilosort output 
% format eg '/mnt/data/N5/2019-09-10_17-01-25'
% output = two structures goodcluster and muacluster for good units and MUA 
% fields:
% (1) id: cluster identity
% (2) spikecount: nb of spikes for this cluster for the full session 
% (3) spikes: spike times in trial 

%% add the relevant repositories to path
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Filepaths
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
addpath (genpath('/mnt/data/N5'))

%session 1 - conc
myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with kilosort output

%session 2 - id
%myKsDir = '/mnt/data/N5/2019-09-11_18-22-53';

%% Get Trial Timestamps from Open Ephys Events file
filename = fullfile(myKsDir,'all_channels.events');
[data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs

% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
offset = offset/sampleRate;
timestamps = timestamps - offset;

%% Get Events 
[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

TrialTimestamps = TrialTimestamps - offset; %adjust for offset
OdorTimestamps = OdorTimestamps - offset; % adjust for offset 

%% Read Stimulus File
%session 1 - conc
stimfilename = '190910_17_01.txt';
%session 2 -id 
%stimfilename = '190911_18_22.txt';

[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList);

%% Load data from kilosort/phy
% fct from spikes package - found in phyhelper

% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unique clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted
sp = loadKSdir(myKsDir);

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
            mygoodspikes = allgoodspikes(find(allgoodspikes>=tstart & allgoodspikes<=tstop));
            mygoodspikes = mygoodspikes - tstart; % align to every trials start
            goodspiketimes(mytrial) = {mygoodspikes}; %spikes in trial
            %each column of goodspiketimes is a trial 
            %for each trial gives you the spiketimes of spikes that occured
        end 
        goodcluster(mycluster).id = sp.cids(mycluster);
        goodcluster(mycluster).spikecount = numel(allgoodspikes);%all spikes
        goodcluster(mycluster).spikes = goodspiketimes; %spikes in trial
        % stores goodspiketimes for each cluster before it gets cleared at
        % next iteration of the loop 
    end
        
    if sp.cgs(mycluster) == 1 % for clusters labeled as MUA
        %get spike times for this cluster 
        allmuaspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds 
        clear muaspiketimes
    
        for mytrial = 1:NTrials
            %TrialSpikes = [];
            tstart = TrialTimestamps(mytrial,1); 
            tstop = TrialTimestamps(mytrial,2); 
            mymuaspikes = allmuaspikes(find(allmuaspikes>=tstart & allmuaspikes<=tstop));
            mymuaspikes = mymuaspikes - tstart; % align every trials 
            muaspiketimes(mytrial) = {mymuaspikes};
        end 
        muacluster(mycluster).id = sp.cids(mycluster);
        muacluster(mycluster).spikecount = numel(allmuaspikes);%all spikes
        muacluster(mycluster).spikes = muaspiketimes; %spikes in trial
    end 
  
end          
disp(['found ',num2str(mycluster),' units']);
disp(['found ',num2str(length (find(sp.cgs==2))),' good units']);
disp(['found ',num2str(length (find(sp.cgs==1))),' multi units']);

end