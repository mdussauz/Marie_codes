%AutocorrelogramCorrelationClosedLoopvsPassive

%% path

% Classic Open Loop:

% Free Lever Sessions:
SessionName = 'S12/2023-08-04_14-29-23';

if strcmp(computer, 'MACI64')
    ephyspath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Ephys/';
else
    ephyspath = '/mnt/data/Processed/Ephys/';
end

myKsDir = fullfile(ephyspath,SessionName);

%% Settings 
to_plot = 1; 
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
    % which tetrode
%     tetrode = floor(sp.channels(mycluster)/4)+1 + ...
%         rem(sp.channels(mycluster),4)/10;

    % Outputs
    % cluster(mycluster).id = sp.cids(mycluster);
    % cluster(mycluster).tetrode = tetrode;
    % cluster(mycluster).spikecount = numel(allspikes);
    % cluster(mycluster).spikes = allspikes;
    % cluster(mycluster).quality = sp.cgs(mycluster);
    % [fpRate, numViolations] = ISIViolations(allspikes, 1/32000, 0.002);
    % cluster(mycluster).ISIquality = [round(fpRate,2,'significant'), round(numViolations/(numel(allspikes)-1),2,'significant')];
%     cluster(mycluster).spikescaling = sp.tempScalingAmps;
%     cluster(mycluster).clusterscalingorder = sp.clu;
%% Extracting spikes and making the ACGs for each part of the session
sessionHalfTime = max(allspikes)/2;

spiketimes_first = allspikes(allspikes <= sessionHalfTime);
spiketimes_second = allspikes(allspikes >= sessionHalfTime);

[xLin1, nLin1, xLog1, nLog1] = myACG(spiketimes_first, [], []);
[xLin2, nLin2, xLog2, nLog2] = myACG(spiketimes_second, [], []);

ACGcorrelation(mycluster) = corr(nLog1/sum(nLog1), nLog2/sum(nLog2),'Type', 'Pearson');

if to_plot 
subnum = mod(mycluster,9);
    if subnum == 0
        subnum = 9;
    end
    subplot(3,3,subnum)
    plot(xLog1,nLog1/(sum(nLog1))); hold on; plot(xLog2, nLog2/(sum(nLog2)))
    title(['unit#',num2str(mycluster),' corr: ',num2str(ACGcorrelation(mycluster))])
    if mod(mycluster,9) == 0
        figure;
    end
end

end

hISI = sgolayfilt(histc(diff(L1chain.times), logspace(-3, 3, 100)), 3, 9) / length(L1chain.times);   