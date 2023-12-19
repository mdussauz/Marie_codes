function [] = PlotWaveformShapeAndDistribution(SessionName, whichunits)
%input:
% SessionName - eg =  'S12/2023-08-04_14-29-23'
%whichunits = units to plot, if not specified, all units will be plotted

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
% sp.cids is list of unique clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted

%% Number of units to plot
if nargin<2
    Nb_Unit = length(sp.cids);
    unit_list = 1:Nb_Unit;
else
    Nb_Unit = length(whichunits);
    unit_list = whichunits;
end

%% get all spikewaveforms
for unit = 1:Nb_Unit % for each cluster
    mycluster =  unit_list(unit);
    % get all spiketimes (in seconds):
    allspikes = sp.st(sp.clu==sp.cids(mycluster));
    % which tetrode:
    tetrode = floor(sp.channels(mycluster)/4)+1 + ...
        rem(sp.channels(mycluster),4)/10;

    % Outputs
    cluster(mycluster).spikescaling = sp.tempScalingAmps; %waveform amplitude
    cluster(mycluster).clusterscalingorder = sp.clu;

    % extract waveforms:
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

    %% Plots
    figure(unit)
    for chunk = 1:3
        whichspikes = intersect(find(wf.spikeTimeKeeps<=1500*chunk),find(wf.spikeTimeKeeps>=1500*(chunk-1)));
        for ch = 1:4
            meanWF = mean(squeeze(wf.waveForms(1,whichspikes,whichchannels(ch),:)),1);
            stdWF = nanstd(squeeze(wf.waveForms(1,whichspikes,whichchannels(ch),:)),1);

            subplot(4,4,ch+(chunk-1)*4);
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
    subplot(4,4,13:16); % amplitude plot
    thisunitamps = cluster(mycluster).spikescaling(find(cluster(mycluster).clusterscalingorder == sp.cids(mycluster)));
    plot(allspikes,thisunitamps,'.');
    set(gca,'YLim',[0 40],'TickDir','out');
    title(['cell# ', num2str(sp.cids(mycluster))])
end
end

        