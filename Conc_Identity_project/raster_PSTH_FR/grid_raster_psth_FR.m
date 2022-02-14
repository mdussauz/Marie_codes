% raster_psth_FR_summary 
% adapted by MD from JW
% Run this code to plot for each cluster of one session: 
% (1) grid of raster plots for each stimulus type and all repeats
% (2) grid of PSTHs averages across repeats for each stimulus type
% (3) grid of FR for odor and air period for each stimulus type
%%
path = '/opt/';
addpath(genpath([path,'Marie_codes']))
addpath(genpath([path,'open-ephys-analysis-tools']));
addpath(genpath([path,'afterphy']));
addpath(genpath([path,'spikes']));
addpath(genpath([path,'npy-matlab']));
addpath(genpath([path, 'Marie_codes/raster PSTH FR']));

dir = '/mnt/data/J2';
experiment = '2020-07-18_16-36-36';
stimfile = '200718_16_36.txt';
concentration = true;

myKsDir = fullfile(dir, experiment);
foo = fullfile(myKsDir, '100_ADC1.continuous');
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs

fast = 0;%1 for dots 0 for lines in raster plots
t_wid = 200; %time window to smooth PSTH over  
% Arka does 100?
%%
% Get offset
OepsSampleRate = 30000; % Open Ephys acquisition rate

[offset] = AdjustClockOffset(myKsDir);
offset = offset/OepsSampleRate;
timestamps = timestamps - offset;

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
        [f] = plotRastersinGrid(cluster,fast); 
        saveas(f, fullfile(dir, 'rasterplots.png'));
        [f] = plotPSTHsinGrid(cluster,t_wid);
        saveas(f, fullfile(dir, 'PSTHs.png'));
        [f] = FiringRateGrid(cluster);
        saveas(f, fullfile(dir, 'firingrate.png'));
        
    end

end