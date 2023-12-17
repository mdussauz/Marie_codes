function [WFcorrelation, conc_WFcorrelation, CorrelationCriteria, concCorrelationCriteria] = GetWaveformCorrelation (SessionName, sessionHalfTime)
%input:
% SessionName - eg =  'S12/2023-08-04_14-29-23'
% sessionHalfTime = Timepoint diving session in two to compare waveform

%% paths
if strcmp(computer, 'MACI64')
    ephyspath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Ephys/';
else
    ephyspath = '/mnt/data/Processed/Ephys/';
end

myKsDir = fullfile(ephyspath,SessionName);

%% Load data from kilosort/phy
sp = loadKSdirPriyanka(myKsDir);
% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unqiue clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted

for mycluster = 1:length(sp.cids) % for each cluster
    % get all spiketimes (in seconds)
    allspikes = sp.st(sp.clu==sp.cids(mycluster));
    % which tetrode
    tetrode = floor(sp.channels(mycluster)/4)+1 + ...
        rem(sp.channels(mycluster),4)/10;

    if nargin<2
        sessionHalfTime = max(allspikes)/2;
    end

    %% get all spikewaveforms
    % get all spikewaveforms
    gwfparams.dataDir = myKsDir;         % KiloSort/Phy output folder
    gwfparams.fileName = sp.dat_path;    % .dat file containing the raw
    gwfparams.dataType = sp.dtype;       % Data type of .dat file (this should be BP filtered)
    gwfparams.nCh = sp.n_channels_dat;   % Number of channels that were streamed to disk in .dat file
    gwfparams.wfWin = [-40 41];          % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = numel(allspikes);    % taking all spikes of a unit; default = 2000; Number of waveforms per unit to pull out - randomly from all spikes of that unit - will take into account if units has less than this number
    gwfparams.spikeTimes = int32(allspikes*sp.sample_rate); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = sp.cids(mycluster) + 0*gwfparams.spikeTimes; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
    wf = getWaveForms(gwfparams);

    wf.spikeTimeKeeps = wf.spikeTimeKeeps/sp.sample_rate;
    whichchannels = (floor(tetrode)*4) + [-3:1:0];

    %% Get Mean Waveform on main channel for each half of session
    mainchannel = rem(sp.channels(mycluster),4)+1;
    %1st half:
    whichspikes_first = find(wf.spikeTimeKeeps<=sessionHalfTime);
    meanWF_first(mycluster,:) = mean(squeeze(wf.waveForms(1,whichspikes_first,whichchannels(mainchannel),:)),1);
    %2nd half
    whichspikes_second = find(wf.spikeTimeKeeps>=sessionHalfTime);
    meanWF_second(mycluster,:) = mean(squeeze(wf.waveForms(1,whichspikes_second,whichchannels(mainchannel),:)),1);

    WFcorrelation(mycluster) = corr(meanWF_first(mycluster,:)', meanWF_second(mycluster,:)','Type', 'Pearson');


    %% Get Concatenated Mean Waveform for th 4 channels
    %init
    conc_meanWF_first = [];
    conc_meanWF_second = [];

    for ch = 1:4
        %1st half:
        whichspikes_first = find(wf.spikeTimeKeeps<=sessionHalfTime);
        temp1 = mean(squeeze(wf.waveForms(1,whichspikes_first,whichchannels(ch),:)),1);
        conc_meanWF_first = horzcat(conc_meanWF_first,temp1);

        %2nd half
        whichspikes_second = find(wf.spikeTimeKeeps>=sessionHalfTime);
        temp2 = mean(squeeze(wf.waveForms(1,whichspikes_second,whichchannels(ch),:)),1);
        conc_meanWF_second = horzcat(conc_meanWF_second, temp2);
    end

    conc_WFcorrelation(mycluster) = corr(conc_meanWF_first', conc_meanWF_second','Type', 'Pearson');

    all_conc_first(mycluster,:) = conc_meanWF_first;


end
    %% Get correlation of the mean waveform on main channel across units for 1st part of session

    pairs = nchoosek(1:length(sp.cids) ,2);
    for comp = 1: length(pairs)
        crossWFcorrelation(comp) = corr(meanWF_first(pairs(comp,1),:)', meanWF_first(pairs(comp,2),:)','Type', 'Pearson');
    end

    CorrelationCriteria = prctile(crossWFcorrelation,99);

    %% Get correlation of the concatenated mean waveforms on the 4 channels across units for first part of session

    pairs = nchoosek(1:length(sp.cids) ,2);
    for comp = 1: length(pairs)
        crossconcWFcorrelation(comp) = corr(all_conc_first(pairs(comp,1),:)', all_conc_first(pairs(comp,2),:)','Type', 'Pearson');
    end

    concCorrelationCriteria = prctile(crossconcWFcorrelation,99);
end
