%raster_psth_FR_summary 

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
foo = fullfile(myKsDir, '100_ADC1.continuous');
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs
%%
% Get offset
OepsSampleRate = 30000; % Open Ephys acquisition rate

[offset] = AdjustClockOffset(myKsDir);
offset = offset/OepsSampleRate;
timestamps = timestamps - offset;
%%
% Get respiration
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