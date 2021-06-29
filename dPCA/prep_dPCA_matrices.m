%prep_dPCA_matrices

%written by MD

% create 2 structure arrays for all neurons all sessions all mice
% 1 structure with spikes aligned to odor onset
% 1 structure with spikes aligned to 1 inhalation after odor onset

%% USER: what brain area

brain_area = 'AON'; % input is 'AON' or 'APC'

%% add the relevant repositories to path
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))
addpath(genpath('/opt/npy-matlab'));
addpath(genpath('/opt/Marie_codes'));
addpath(genpath('/opt/Marie_codes/dPCA'))

%% Filepaths
%folder with stimulus files
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
%folders with sessions coming from AON mice
addpath (genpath('/mnt/interim/N01'))
addpath (genpath('/mnt/interim/N02'))
addpath (genpath('/mnt/interim/D1'))
addpath (genpath('/mnt/interim/D2'))
addpath (genpath('/mnt/interim/D3'))
%folders with sessions coming from APC mice
addpath (genpath('/mnt/interim/PJ1'))
addpath (genpath('/mnt/interim/PJ2'))
addpath (genpath('/mnt/interim/APC1'))
addpath (genpath('/mnt/interim/APC2'))

%% Recording sessions

if brain_area == 'AON'
    myKsDir = {'/mnt/interim/N01/2021-01-25_09-35-25','/mnt/interim/N01/2021-01-27_10-46-51',...
        '/mnt/interim/N02/2021-01-26_11-47-01','/mnt/interim/N02/2021-01-28_10-57-31',...
        '/mnt/interim/N02/2021-01-29_14-06-12',...
        '/mnt/interim/D1/2021-04-27_09-30-31', '/mnt/interim/D1/2021-05-05_11-42-47',...
        '/mnt/interim/D1/2021-05-06_13-53-52', '/mnt/interim/D1/2021-05-07_09-13-39',...
        '/mnt/interim/D3/2021-04-26_17-24-41', '/mnt/interim/D3/2021-04-28_15-59-15',...
        '/mnt/interim/D3/2021-04-29_16-43-06','/mnt/interim/D3/2021-04-30_16-23-08',...
        '/mnt/interim/D3/2021-05-20_13-54-53','/mnt/interim/D3/2021-05-21_13-04-28',...
        '/mnt/interim/D3/2021-05-24_14-33-50','/mnt/interim/D3/2021-05-25_14-39-17',...
        '/mnt/interim/D3/2021-05-26_13-18-29','/mnt/interim/D3/2021-05-27_14-11-45',...
        '/mnt/interim/D2/2021-05-31_13-19-08', '/mnt/interim/D2/2021-06-01_12-29-35',...
        '/mnt/interim/D2/2021-06-02_13-28-26', '/mnt/interim/D2/2021-06-03_15-09-14',...
        '/mnt/interim/D2/2021-06-09_12-39-32','/mnt/interim/D2/2021-06-10_12-56-18'...
        '/mnt/interim/D2/2021-06-15_15-43-57', '/mnt/interim/D2/2021-06-16_10-35-59'};
elseif brain_area == 'APC'
    myKsDir = {'/mnt/interim/PJ1/2021-01-25_15-13-27','/mnt/interim/PJ1/2021-01-28_14-21-40',...
        '/mnt/interim/PJ2/2021-02-14_15-06-20','/mnt/interim/PJ2/2021-02-15_14-51-10',...
        '/mnt/interim/APC1/2021-04-26_11-12-59', '/mnt/interim/APC1/2021-04-27_15-37-37',...
        '/mnt/interim/APC1/2021-04-28_09-23-05','/mnt/interim/APC1/2021-04-29_10-12-26',...
        '/mnt/interim/APC1/2021-05-03_14-58-49', '/mnt/interim/APC1/2021-05-04_14-16-04',...
        '/mnt/interim/APC2/2021-04-27_12-35-49', '/mnt/interim/APC2/2021-04-28_12-09-32',...
        '/mnt/interim/APC2/2021-04-29_13-19-18', '/mnt/interim/APC2/2021-04-30_13-26-00',...
        '/mnt/interim/APC2/2021-05-03_11-12-06', '/mnt/interim/APC2/2021-05-04_11-24-28'...
        };
end

%% Stimulus File

if brain_area == 'AON'
    stimfilename = {'210125_9_35.txt', '210127_10_47.txt',... %N01
        '210126_11_47.txt','210128_10_58.txt', '210129_14_06.txt',...%N02
        '210427_9_30.txt','21055_11_42.txt','21056_13_53.txt','21057_9_13.txt',... %D1
        '210426_17_24.txt','210428_15_58.txt','210429_16_42.txt','210430_16_22.txt',...%D3
        '210520_13_55.txt', '210521_13_04.txt', '210524_14_33.txt',...
        '210525_14_39.txt','210526_13_18.txt', '210527_14_11.txt',...
        '210531_13_19.txt', '21061_12_29.txt', '21062_13_33.txt', '21063_15_09.txt', ...%D2 %'21062_13_28.txt' is incorrect
        '21069_12_39.txt','210610_12_56.txt','210615_15_44.txt','210616_10_36.txt'...
        };
elseif brain_area == 'APC'
    stimfilename = {'210125_15_13.txt', '210128_14_21.txt',... %PJ1
        '210214_15_06.txt', '210215_14_52.txt',... %PJ2
        '210426_11_12.txt', '210427_15_37.txt', '210428_9_22.txt',...%APC1
        '210429_10_11.txt','210430_9_23.txt','21053_14_58.txt','21054_14_15.txt',...
        '210427_12_35.txt', '210428_12_09.txt','210429_13_18.txt', ... %APC2
        '210430_13_25.txt','21053_11_11.txt','21054_11_23.txt',...
        };
end

%% defaults and initialize structures 
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz

%if exist('allclusters','var')== 0 && exist('allrespiclusters','var')==0
allclusters = struct('id',[],'spikecount',[],'spikes',[], 'stimulus', []); %initiate
allrespiclusters = struct('id',[],'spikecount',[],'spikes',[], 'stimulus', []); %initiate
%end 

%% Get the structures  

for i = 1: length(myKsDir) % loop through each recording session
    %% Get Trial Timestamps from Open Ephys Events file
    [StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename{i}); %get stimulus list from stimulus file
    NOdors = numel(unique(StimList));
    NTrials= numel(StimList);
    disp(['found ',num2str(NTrials),' trials']); % for checking purpose
    
    filename = fullfile(myKsDir{i},'all_channels.events');
    [data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs
    
    % adjust for clock offset between open ephys and kilosort
    [offset] = AdjustClockOffset(myKsDir{i});
    offset = offset/sampleRate;
       
    %% Get timestamps of Events (trial on/off and odor on/off) and correct for ephys offset
    [Events] = ParseOpenEphysEvents(filename); 
    % trial events are in Channel0, odor events are in Channel1
    
   % if filename == '/mnt/interim/PJ2/2021-02-14_15-06-20/all_channels.events'
   if strcmpi('/mnt/interim/PJ2/2021-02-14_15-06-20/all_channels.events',filename)
        Events.Channel0.Off(1) = [];
    end
    
    TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
    TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
    OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
    
    TrialTimestamps = TrialTimestamps - offset; %adjust for offset
    OdorTimestamps = OdorTimestamps - offset; % adjust for offset
    
    %% Get respiration onset timestamps and correct for ephys offset
    [closest_neighbor_zci] = getFirstInhalation(myKsDir{i}, OdorTimestamps);
    
    %% Load data from kilosort/phy
    sp = loadKSdir(myKsDir{i}); % get all units
    
    %% Number of good units - for info 
    
    goodUNber = length (find(sp.cgs==2)); % count how many good units
    disp(['found ',num2str(goodUNber),' good units']); % display for user
    
    %% Split data by clusters and by trials
    
    for mycluster = 1:length(sp.cids) % for each cluster
        
        if sp.cgs(mycluster) == 2 % for clusters labeled as good
            %get spike times for this cluster
            allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds
            clear goodspiketimes
            clear respspiketimes
            
            for mytrial = 1:NTrials % for each trial 
            
                tstart = TrialTimestamps(mytrial,1); %start of trial
                tstop = TrialTimestamps(mytrial,2); %end of trial
                odorstart = OdorTimestamps(mytrial,1); % odor onset 
                inhalationStart = closest_neighbor_zci(mytrial); % 1st inhalation after odor onset
                
                mygoodspikes = allgoodspikes(find(allgoodspikes>=tstart & allgoodspikes<=tstop));% only keep spikes that are in trial
                mygoodspikes = mygoodspikes - odorstart; % align to every odor start
                myrespirationspikes = mygoodspikes - inhalationStart; % align to every 1st inhalation
                
                goodspiketimes(mytrial) = {mygoodspikes}; %timestamps of spikes in trial and aligned to odor start
                respspiketimes(mytrial) = {myrespirationspikes}; %timestamps of spikes in trial and aligned to 1st inhalation
                %each column of goodspiketimes is a trial
                %for each trial gives you the spiketimes of spikes that occured
            end
            
            goodcluster(mycluster).id = sp.cids(mycluster); % cluster identity
            goodcluster(mycluster).spikecount = numel(allgoodspikes);%count of all spikes in trial
            goodcluster(mycluster).spikes = goodspiketimes; %timestamps of spikes in trial
            % stores goodspiketimes for each cluster before it gets cleared at
            % next iteration of the loop
            goodcluster(mycluster).stimulus = StimList; %list of stimuli
            
            respicluster(mycluster).id = sp.cids(mycluster);
            respicluster(mycluster).spikecount = numel(allgoodspikes);%all spikes
            respicluster(mycluster).spikes = respspiketimes; %spikes in trial
            % stores respspiketimes for each cluster before it gets cleared at
            % next iteration of the loop
            respicluster(mycluster).stimulus = StimList; %list of stimuli 
        end
        
    end
    
    allclusters =  [allclusters, goodcluster];
    allclusters( all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), allclusters, 'UniformOutput', false ) ), 1 ) ) = []; %remove empty lines
    
    allrespiclusters = [allrespiclusters,respicluster];
    allrespiclusters( all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), allrespiclusters, 'UniformOutput', false ) ), 1 ) ) = []; %remove empty lines
    
    clear goodcluster ; clear respicluster
end

%save ('/opt/allclusters.m',  'allclusters')