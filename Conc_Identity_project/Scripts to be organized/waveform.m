%% Filepaths
addpath (genpath('/mnt/data/D4'))

%session 1 - conc
myKsDir = '/mnt/data/D4/2021-07-23_10-27-41'; % directory with kilosort output

%% Load data from kilosort/phy
% fct from spikes package - found in phyhelper

% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unique clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted
sp = loadKSdir(myKsDir);

%%
for mycluster = 1:length(sp.cids) % for each cluster 
    
    if sp.cgs(mycluster) == 2 % for clusters labeled as good
        %get spike times for this cluster 
        allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds

%% INPUT

gwfparams.dataDir = myKsDir;             % KiloSort/Phy output folder
gwfparams.fileName = 'mybinaryfile.dat'; % .dat file containing the raw 
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

%%
wf = getWaveForms(gwfparams);

%%
