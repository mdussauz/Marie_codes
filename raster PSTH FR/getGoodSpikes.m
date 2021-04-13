function [goodcluster, StimList, TrialTimestamps, OdorTimestamps] = getGoodSpikes(myKsDir, offset, stimfile)
% gets good clusters
% Input: 
% - myKsDir = directory with data
% - offset =  offset
% stim_file = the name of the stimulus file, ex: "191021.txt"

%%
path = '~/opt/';
addpath(genpath([path,'open-ephys-analysis-tools']));
addpath(genpath([path,'afterphy']));
addpath(genpath([path,'spikes']));
addpath(genpath([path,'npy-matlab']));

%%
filename = fullfile(myKsDir,'all_channels.events');
[data, timestamps, info] = load_open_ephys_data(filename);
[Events] = ParseOpenEphysEvents(filename);

stimfilename = fullfile(myKsDir, stimfile);
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList)/NOdors;
% trial events are in Channel0, odor events are in Channel1

%%
TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

TrialTimestamps = TrialTimestamps - offset; %adjust for offset
OdorTimestamps = OdorTimestamps - offset;

%%
sp = loadKSdir(myKsDir);

for mycluster = 1:length(sp.cids) % for each cluster 
    
    if sp.cgs(mycluster) == 2 % for clusters labeled as good
        %get spike times for this cluster 
        allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds 
        
        goodspiketimes = cell(NOdors, NTrials);
        odorspiketimes = cell(NOdors, NTrials);
        trial = 1;
        for ind = 1:length(StimList)
            %TrialSpikes = [];
            tstart = TrialTimestamps(ind,1); 
            tstop = TrialTimestamps(ind,2);  
            odorstart = OdorTimestamps(ind,1); 
            mygoodspikes = allgoodspikes(find(allgoodspikes>=tstart & allgoodspikes<=tstop));
            odorgoodspikes = mygoodspikes - odorstart; % align to every odor start!!!
            goodspiketimes(StimList(ind), trial) = {mygoodspikes};
            odorspiketimes(StimList(ind), trial) = {odorgoodspikes};%spikes in trial
            if mod(ind, NOdors) == 0
                trial = trial + 1;
            end
            %each column of goodspiketimes is a trial 
            %for each trial gives you the spiketimes of spikes that occured
        end 
        goodcluster(mycluster).id = sp.cids(mycluster);
        goodcluster(mycluster).spikecount = numel(allgoodspikes);%all spikes
        goodcluster(mycluster).spikes = goodspiketimes;%spikes in trial
        goodcluster(mycluster).odorspikes = odorspiketimes; %spikes aligned to odor start
        % stores goodspiketimes for each cluster before it gets cleared at
        % next iteration of the loop 
    end
        
  
end



end
