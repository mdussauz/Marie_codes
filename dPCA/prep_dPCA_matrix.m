%dPCA_matrix

%written by MD

% create structure array for all neurons all sessions all mice


%% add the relevant repositories to path
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))
addpath(genpath('/opt/npy-matlab'));
addpath(genpath('/opt/Marie_codes'));
addpath(genpath('/opt/Marie_codes/dPCA'))

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Stimulus File
% stimfilename = {'190910_17_01.txt', '190918_16_35.txt', '190921_17_29.txt', ...
%     '190925_13_28.txt', '190912_14_01.txt', '190919_13_40.txt', ...
%     '190923_13_48.txt', '190925_11_39.txt', '190910_15_23.txt', ...
%     '190912_15_42.txt', '190913_15_47.txt', '190919_12_05.txt', ...
%     '190923_17_30.txt','200718_16_36.txt'};
stimfilename = {'210125_9_35.txt', '210127_10_47.txt', '210126_11_47.txt', ...
     '210128_10_58.txt', '210129_14_06.txt'};

%% Filepaths
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
% addpath (genpath('/mnt/data/N5'))
% addpath (genpath('/mnt/data/N6'))
% addpath (genpath('/mnt/data/N8'))
% addpath (genpath('/mnt/data/N10'))
% addpath (genpath('/mnt/data/J2'))
addpath (genpath('/mnt/data/N01'))
addpath (genpath('/mnt/data/N02'))

% myKsDir = {'/mnt/data/N5/2019-09-10_17-01-25','/mnt/data/N5/2019-09-18_16-35-17',...
%     '/mnt/data/N5/2019-09-21_17-29-13', '/mnt/data/N5/2019-09-25_13-28-14', ...
%     '/mnt/data/N6/2019-09-12_14-01-04', '/mnt/data/N6/2019-09-19_13-40-40', ...
%     '/mnt/data/N6/2019-09-23_13-48-31', '/mnt/data/N6/2019-09-25_11-39-08', ...
%     '/mnt/data/N8/2019-09-10_15-22-54', '/mnt/data/N8/2019-09-12_15-42-26', ...
%     '/mnt/data/N10/2019-09-13_15-47-27', '/mnt/data/N10/2019-09-19_12-04-11', ...
%     '/mnt/data/N10/2019-09-23_17-30-17','/mnt/data/J2/2020-07-18_16-36-36'};

myKsDir = {'/mnt/data/N01/2021-01-25_09-35-25','/mnt/data/N01/2021-01-27_10-46-51',...
    '/mnt/data/N02/2021-01-26_11-47-01','/mnt/data/N02/2021-01-28_10-57-31',...
    '/mnt/data/N02/2021-01-29_14-06-12'};


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

%save ('/opt/allclusters.m',  'allclusters')