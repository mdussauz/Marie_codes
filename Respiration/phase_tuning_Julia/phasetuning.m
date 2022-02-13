
path = '~/Desktop/code/';
addpath(genpath([path,'open-ephys-analysis-tools']));
addpath(genpath([path,'afterphy']));
addpath(genpath([path,'spikes']));
addpath(genpath([path,'npy-matlab']));
addpath(genpath([path, 'phase tuning']));

dir = '~/Desktop/Analysis/';
experiment = '2019-09-10_17-01-25';
stimfile = '190910_17_01.txt';
concentration = true;

myKsDir = fullfile(dir, experiment);
foo = fullfile(myKsDir, '100_ADC1.continuous'); % pressure sensor
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs
%%
% Get offset
% adjust timestamps to account for the start offset in OEPS
OepsSampleRate = 30000; % Open Ephys acquisition rate

[offset] = AdjustClockOffset(myKsDir);
offset = offset/OepsSampleRate;
timestamps = timestamps - offset;
%%
% Get respiration
% downsample respiration data to 1 KHz - necessary
newTimeBin = 0.001; % in seconds
newTimeBase = (0:newTimeBin:max(timestamps))';
Respiration_DS = interp1q(timestamps,Respiration,newTimeBase);
%%
% plot respiration and get respiration phase starts
[MyPeaks, Idx, Respiration_DS_Filt] = GetRespirationTimeStamps(Respiration_DS,'plotfigures',1,'threshold',0.04);
phaseStarts = newTimeBase(Idx);

%%
% get good clusters
[goodcluster, StimList] = getGoodSpikes(myKsDir, offset, stimfile);
%%
% make graphs
for cluster = goodcluster
    if ~isempty(cluster.id)
        dir = fullfile(myKsDir, strcat("cluster", num2str(cluster.id)));
        if exist(dir, 'dir')
            rmdir(dir, 's')
        end
        mkdir(dir);
        f = plotSpikesinPhaseGrid(cluster, phaseStarts);
        saveas(f, fullfile(dir, 'phaseraster.png'));
        [f, aircurves, odorcurves] = phaseCurvesGrid(cluster, phaseStarts, 10);
        saveas(f, fullfile(dir, 'tuningcurves.png'));
        f = correlationHeatMap(aircurves, odorcurves);
        saveas(f, fullfile(dir, 'corrheatmap.png'));
        [f, deltaThetas] = plotPolarGrid(cluster, phaseStarts);
        saveas(f, fullfile(dir, 'polarplots.png'));
        [f] = FiringRateGrid(cluster);
        saveas(f, fullfile(dir, 'firingrate.png'));
    end

end

%%
% rest is for summary figure
mean_air_theta = [];
mean_air_rho = [];
mean_odor_theta = [];
mean_odor_rho = [];
%%
for cluster = goodcluster
    if ~isempty(cluster.id)
        [airthetas, airrhos, odorthetas, odorrhos] = getAllVectors(cluster, phaseStarts);
        mean_air_theta = vertcat(mean_air_theta, mean(airthetas, 1));
        mean_air_rho = vertcat(mean_air_rho, mean(airrhos, 1));
        mean_odor_theta = vertcat(mean_odor_theta, mean(odorthetas, 1));
        mean_odor_rho = vertcat(mean_odor_rho, mean(odorrhos, 1));
    end
end
%%
figure()
subplot(1, 2, 1)
for cluster = 1:size(mean_air_rho, 1)
   polarplot([0 mean_air_theta(cluster, 1)], [0 mean_air_rho(cluster, 1)])
   hold on
end
title("air")
rlim([0 1])
subplot(1, 2, 2)
for cluster = 1:size(mean_odor_rho, 1)
   polarplot([0 mean_odor_theta(cluster, 1)], [0 mean_odor_rho(cluster, 1)])
   hold on
end
rlim([0 1])
title("odor")
sgtitle("Concentration 10-4");