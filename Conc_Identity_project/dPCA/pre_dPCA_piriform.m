%dPCA_matrix

%written by MD

% create structure array for all neurons all sessions all mice


%% add the relevant repositories to path
addpath(genpath('/opt/Marie_codes'))
addpath(genpath('/opt/Marie_codes/dPCA'))
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))
addpath(genpath('/opt/npy-matlab'));

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Stimulus File
% stimfilename = {'200623_15_17.txt', '200625_15_14.txt', '200721_14_25.txt', ...
%     '20087_16_33.txt', '200713_11_26.txt', '200718_18_34.txt', ...
%     '20087_14_48.txt', '200810_12_24.txt'};

stimfilename = {'210125_15_13.txt', '210128_14_21.txt', '210214_15_06.txt', ...
     '210215_14_52.txt'};
%% Filepaths
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
% addpath (genpath('/mnt/data/J3'))
% addpath (genpath('/mnt/data/J4'))
addpath (genpath('/mnt/data/PJ1'))
addpath (genpath('/mnt/data/PJ2'))

% myKsDir = {'/mnt/data/J3/2020-06-23_15-17-20','/mnt/data/J3/2020-06-25_15-13-48',...
%     '/mnt/data/J3/2020-07-21_14-25-37', '/mnt/data/J3/2020-08-07_16-33-05', ...
%     '/mnt/data/J4/2020-07-13_11-26-26', '/mnt/data/J4/2020-07-18_18-33-57', ...
%     '/mnt/data/J4/2020-08-07_14-48-11', '/mnt/data/J4/2020-08-10_12-24-03'};

myKsDir = {'/mnt/data/PJ1/2021-01-25_15-13-27','/mnt/data/PJ1/2021-01-28_14-21-40',...
    '/mnt/data/PJ2/2021-02-14_15-06-20','/mnt/data/PJ2/2021-02-15_14-51-10'};


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

if filename == '/mnt/data/PJ2/2021-02-14_15-06-20/all_channels.events'
    Events.Channel0.Off(1) = [];
end
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