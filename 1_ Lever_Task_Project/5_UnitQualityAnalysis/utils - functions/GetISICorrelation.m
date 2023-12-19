function [ISIcorrelation] = GetISICorrelation(SessionName, sessionHalfTime)
%input:
% SessionName - eg =  'S12/2023-08-04_14-29-23'
% sessionHalfTime = Timepoint diving session in two to compare waveform

%% path
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

%%
for mycluster = 1:length(sp.cids) % for each cluster
    % get all spiketimes (in seconds)
    allspikes = sp.st(sp.clu==sp.cids(mycluster));

    if nargin<2
        sessionHalfTime = max(allspikes)/2;
    end

    %% Extracting spikes and making the ACGs for each part of the session

    spiketimes_first = allspikes(allspikes <= sessionHalfTime);
    spiketimes_second = allspikes(allspikes >= sessionHalfTime);

    hISI1 = sgolayfilt(histc(diff(spiketimes_first), logspace(-3, 3, 100)), 3, 9) / length(spiketimes_first);
    hISI2 = sgolayfilt(histc(diff(spiketimes_second), logspace(-3, 3, 100)), 3, 9) / length(spiketimes_second);

    ISIcorrelation(mycluster) = corr(hISI1, hISI2,'Type', 'Pearson');


end
end