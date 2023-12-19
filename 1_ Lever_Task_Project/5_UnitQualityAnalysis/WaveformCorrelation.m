%WaveformCorrelation

% From Schoonover et al 2021
% Single-unit waveform similarity is defined as the Pearson’s correlation
% between a pair of concatenated average waveform vectors.
% Within a given recording day, the 99th percentile correlation between
% the waveforms of different single units was 0.93.
% We therefore excluded single units whose correlation across a pair of
% comparison days fell below this value. Variations in waveform shape are
% amplified in large-amplitude single units. Thus, we used a final manual
% check to identify and include those rare cases in which very large
% amplitude single units had across-day waveform correlations that
% fell slightly below 0.93

%% Settings
to_plot = 0; % 0 to skip plotting
%% paths

% Classic Open Loop:

% Free Lever Sessions:
SessionName = 'S12/2023-08-04_14-29-23';

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
    % which tetrode
    tetrode = floor(sp.channels(mycluster)/4)+1 + ...
        rem(sp.channels(mycluster),4)/10;

    % Outputs
    % cluster(mycluster).id = sp.cids(mycluster);
    % cluster(mycluster).tetrode = tetrode;
    % cluster(mycluster).spikecount = numel(allspikes);
    % cluster(mycluster).spikes = allspikes;
    % cluster(mycluster).quality = sp.cgs(mycluster);
    % [fpRate, numViolations] = ISIViolations(allspikes, 1/32000, 0.002);
    % cluster(mycluster).ISIquality = [round(fpRate,2,'significant'), round(numViolations/(numel(allspikes)-1),2,'significant')];
    cluster(mycluster).spikescaling = sp.tempScalingAmps;
    cluster(mycluster).clusterscalingorder = sp.clu;
    
    %% get all spikewaveforms
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
    sessionHalfTime = max(allspikes)/2;
    %sessionHalfTime = 4000;
    mainchannel = rem(sp.channels(mycluster),4)+1;
    %1st half:
    whichspikes_first = find(wf.spikeTimeKeeps<=sessionHalfTime);
    meanWF_first(mycluster,:) = mean(squeeze(wf.waveForms(1,whichspikes_first,whichchannels(mainchannel),:)),1);
    %2nd half
    whichspikes_second = find(wf.spikeTimeKeeps>=sessionHalfTime);
    meanWF_second(mycluster,:) = mean(squeeze(wf.waveForms(1,whichspikes_second,whichchannels(mainchannel),:)),1);

    WFcorrelation(mycluster) = corr(meanWF_first(mycluster,:)', meanWF_second(mycluster,:)','Type', 'Pearson');
    
    subnum = mod(mycluster,9);
    if subnum == 0
        subnum = 9;
    end
    subplot(3,3,subnum)
    plot(meanWF_first(mycluster,:)); hold on;plot(meanWF_second(mycluster,:))
    title(['unit#',num2str(mycluster),' corr: ',num2str(WFcorrelation(mycluster))])
    if mod(mycluster,9) == 0
        figure;
    end
    
    %% Get Concatenated Mean Waveform for th 4 channels

    sessionHalfTime = 4000;
    mainchannel = rem(sp.channels(mycluster),4)+1;

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
    conc_first(mycluster,:) = conc_meanWF_first;
    conc_second(mycluster,:) = conc_meanWF_second;

    conc_WFcorrelation(mycluster) = corr(conc_meanWF_first', conc_meanWF_second','Type', 'Pearson');


    %% Get correlation across units

    pairs = nchoosek(1:length(sp.cids) ,2);
    for comp = 1: length(pairs)
        crossWFcorrelation(comp) = corr(meanWF_first(pairs(comp,1),:)', meanWF_first(pairs(comp,2),:)','Type', 'Pearson');
    end
    %% Plotting
    if to_plot
        figure;
        for chunk = 1:3
            whichspikes = intersect(find(wf.spikeTimeKeeps<=1500*chunk),find(wf.spikeTimeKeeps>=1500*(chunk-1)));
            for ch = 1:4
                subplot(4,4,ch+(chunk-1)*4);
                meanWF = mean(squeeze(wf.waveForms(1,whichspikes,whichchannels(ch),:)),1);
                stdWF = nanstd(squeeze(wf.waveForms(1,whichspikes,whichchannels(ch),:)),1);


                plot(meanWF,'k');
                hold on
                plot(meanWF+stdWF,'b');
                plot(meanWF-stdWF,'b');
                semWF = stdWF/sqrt(numel(whichspikes));
                plot(meanWF+semWF,'r');
                plot(meanWF-semWF,'r');
                %MyShadedErrorBar(1:length(meanWF),meanWF,stdWF,'r');
                set(gca,'YLim',[-600 250],'XLim',[0 83],'TickDir','out');

            end
        end
    end
    if to_plot
        subplot(4,4,[13:16]); % amplitude plot
        thisunitamps = cluster(mycluster).spikescaling(find(cluster(mycluster).clusterscalingorder == sp.cids(mycluster)));
        plot(allspikes,thisunitamps,'.');
        set(gca,'YLim',[0 40],'TickDir','out');
    end
end

