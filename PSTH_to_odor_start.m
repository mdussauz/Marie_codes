%PSTH_to_odor_start 

%written by MD

%mycluster = 1;

clear all 

%% add the relevant repositories to path
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Read Stimulus File

% MOUSE N5 
%session 1 - conc
%stimfilename = '190910_17_01.txt';
%session 2 -id 
%stimfilename = '190911_18_22.txt';
%session 2 - conc 
stimfilename = '190918_16_35.txt';
%session 3 - conc 
%stimfilename = '190921_17_29.txt';
%session 4 - conc 
%stimfilename = '190925_13_28.txt';

% MOUSE N6 
%session 1 - conc 
%stimfilename = '190909_14_02.txt'; % doesn't work 
%session 2 - conc 
%stimfilename = '190912_14_01.txt';
%session 3 - conc 
%stimfilename = '190919_13_40.txt';
%session 4 - conc 
%stimfilename = '190921_15_38.txt'; doesn't work 
%session 5 - conc 
%stimfilename = '190923_13_48.txt'; 
%session 6 - conc 
%stimfilename = '190925_11_39.txt'; 

% MOUSE N8
%session 1 - conc
%stimfilename = '190910_15_23.txt';
%session 2 - conc
%stimfilename = '190912_15_42.txt';

% MOUSE N10
%session 1 - conc
%stimfilename = '190913_15_47.txt';
%session 2 - conc
%stimfilename = '190919_12_05.txt';
%session 3 - conc
%stimfilename = '190923_17_30.txt';

[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList);

%% Filepaths
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))

% MOUSE N5 
addpath (genpath('/mnt/data/N5'))

%session 1 - conc
%myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with kilosort output
%session 2 - id
%myKsDir = '/mnt/data/N5/2019-09-11_18-22-53';
%session 2 - conc
myKsDir = '/mnt/data/N5/2019-09-18_16-35-17'; % directory with kilosort output
%session 3 - conc
%myKsDir = '/mnt/data/N5/2019-09-21_17-29-13'; % directory with kilosort output
%session 4 - conc
%myKsDir = '/mnt/data/N5/2019-09-25_13-28-14'; % directory with kilosort output

% MOUSE N6 
%addpath (genpath('/mnt/data/N6'))

%session 1 - conc
%myKsDir = '/mnt/data/N6/2019-09-09_14-02-33'; % doesn't work 
%session 2 - conc
%myKsDir = '/mnt/data/N6/2019-09-12_14-01-04'; % directory with kilosort output
%session 3 - conc
%myKsDir = '/mnt/data/N6/2019-09-19_13-40-40'; % directory with kilosort output
%session 4 - conc
%myKsDir = '/mnt/data/N6/2019-09-21_15-38-36'; % doesn't work 
%session 5 - conc
%myKsDir = '/mnt/data/N6/2019-09-23_13-48-31'; % directory with kilosort output
%session 6 - conc
%myKsDir = '/mnt/data/N6/2019-09-25_11-39-08'; % directory with kilosort output

% MOUSE N8
%addpath (genpath('/mnt/data/N8'))

%session 1 - conc
%myKsDir = '/mnt/data/N8/2019-09-10_15-22-54'; % directory with kilosort output
%session 2 - conc
%myKsDir = '/mnt/data/N8/2019-09-12_15-42-26'; % directory with kilosort output

% MOUSE N10
%addpath (genpath('/mnt/data/N10'))

%session 1 - conc
%myKsDir = '/mnt/data/N10/2019-09-13_15-47-27'; % directory with kilosort output
%session 2 - conc
%myKsDir = '/mnt/data/N10/2019-09-19_12-04-11'; % directory with kilosort output
%session 3 - conc
%myKsDir = '/mnt/data/N10/2019-09-23_17-30-17'; % directory with kilosort output

%% Get Trial Timestamps from Open Ephys Events file
filename = fullfile(myKsDir,'all_channels.events');
[data, timestamps, info] = load_open_ephys_data(filename); % data has channel IDs

% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
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
% fct from spikes package - found in phyhelper

% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unique clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted
sp = loadKSdir(myKsDir);

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
    end
        
  
end          

%%
for whichcluster = 1:length(goodcluster)
    clear spikecount
  if ~isempty(goodcluster(whichcluster).id)    
% PSTH of one unit for all trials %%%%

step = 0.001; %timebin for PSTH should be 1 ms 

for mytrial = 1:NTrials
%    mytrial 
%     odorstart = OdorTimestamps(mytrial,1);
% 	tstart = TrialTimestamps(mytrial,1) - odorstart; 
% 	tstop = TrialTimestamps(mytrial,2) - odorstart; 
%   lastbin = ceil(tstop)-1
    firstbin = -10;
    lastbin = 10;
    tbins = firstbin:step:lastbin;
    SpikeTimes = goodcluster(whichcluster).spikes{1, mytrial};
    myspikecount = histcounts(SpikeTimes, tbins);
    spikecount(mytrial,:) = [myspikecount];
end

meanalltrialspikecount = ((sum(spikecount, 1)).*(1/step))./NTrials; %sum(A,1) = sum along each column of the matrix A 

taxis = -500:500;  % make a time axis of 1000 ms
t_wid = 100;  % width of kernel (in ms)
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

%clusterNum = 1:29;

% for c = 1:length(clusterNum) % loop through each unit
%     clusterIdx = clusterNum(c);
%     for t = 1:size(NTrials,1) % loop through trials
%         tempPSTH = squeeze(PSTH(whichcluster,t,:));

        %PSTH is where I store spike times (1 for spikes and 0 elsewhere)
        %PSTH is Neurons X trials X time (in ms)

        zs = conv(meanalltrialspikecount,gauss_kernel,'same');
        FR = zs*1000; %converting firing rate to Hz.

        % FR contains smoothed firing rates;

%    end
% end

smoothmeanall = smoothdata(meanalltrialspikecount, 'gaussian', 10);
% %test = mean(spikecount); %gets mean of each column 
% stdalltrialspikecount = std((sum(spikecount, 1)).*(1/step)) ;
% semalltrialspikecount = stdalltrialspikecount ./sqrt(NTrials);

figure(whichcluster)
plot(tbins(1:end-1),meanalltrialspikecount, '-k');  hold on
plot(tbins(1:end-1),smoothmeanall, '-c', 'LineWidth', 2);  hold on

%xlim([0 max(tbins)]);
xlim([-10 10]);
%plot([10 10], [0 max(meanalltrialspikecount)+1], 'r--'); hold on % dotted red line odor on
%plot([14 14], [0 max(meanalltrialspikecount)+1], 'r--'); hold on
title(['PSTH of cluster #', num2str(goodcluster(whichcluster).id), ' for all trials', ' - tbins(s) ',  num2str(step) ])
hold off

figure(whichcluster + goodUNber)
plot(tbins(1:end-1), zs, '-c', 'LineWidth', 2);


%
%figure(3); imagesc(meanalltrialspikecount);

% %% PSTH of one unit sorted by concentation and odor identity 
% 
% figure(whichcluster + goodUNber)
% sgtitle(['PSTH of cluster #', num2str(goodcluster(whichcluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step) ])
% 
% for i = 1:NOdors %for each odor
%     reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
%     clear spikecount 
% 
%     
%     for j = 1:numel(reps)
%             
%     tstart = TrialTimestamps(reps(j),1); 
% 	tstop = TrialTimestamps(reps(j),2);  
%     lastbin = ceil(tstop-tstart); 
%     tbins = 0:step:lastbin; 
%     SpikeTimes = goodcluster(whichcluster).spikes{1, reps(j)};
%     myspikecount = histcounts(SpikeTimes, tbins);
%     spikecount(j,:) = [myspikecount];
%     
% 
%     end
%         
%         meanspikecount = ((sum(spikecount, 1)).*(1/step))./numel(reps);
%         smoothmean = smoothdata(meanspikecount, 'gaussian', 10);
%         subplot(4, 5, StimList(i))
%         plot(tbins(1:end-1),meanspikecount, '-k', 'LineWidth', 0.4); hold on
%         plot(tbins(1:end-1),smoothmean, '-c', 'LineWidth', 1);  hold on
% %         figure(5); subplot(4, 5, StimList(i));  imagesc(meanspikecount);
%         
%         clear meanspikecount
%  end 
%     
% %totalspike = sum(allspikecount,1)
%     for k = 1:20
%         
%         subplot(4, 5, k)
%         
%        % plot([10 10], [0 40], 'r--'); hold on % dotted red line odor on
%        % plot([14 14], [0 40], 'r--'); hold on % dotted red line odor off
%         
%         %axis
%         %xlim([-20 21])
%         xlim([0 10]);
%         ylim([0 max(meanalltrialspikecount)+1])
%        
%     end
%     

  end
end
