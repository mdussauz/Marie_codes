

path = '~/Desktop/code/';
addpath(genpath([path,'Marie']))
addpath(genpath([path,'open-ephys-analysis-tools']));
addpath(genpath([path,'afterphy']));
addpath(genpath([path,'spikes']));
addpath(genpath([path,'npy-matlab']));
addpath(genpath([path, 'channel and ADC check']));

dir = '~/Desktop/Analysis/';
experiment = '2019-09-10_17-01-25';
stimfile = '190910_17_01.txt';
concentration = true;
channelNum = 32; %nb of channels - change to 64 for 16 TT

myKsDir = fullfile(dir, experiment);

%% Respiration data 
% extract respiration data 
foo = fullfile(myKsDir, '100_ADC1.continuous');
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs

% adjust timestamps to account for the start offset in OEPS
OepsSampleRate = 30000; % Open Ephys acquisition rate
%[offset] = AdjustClockOffset(myKsDir);
offset = timestamps(1);
timestamps = timestamps - offset;

% downsample respiration data to 1 KHz - necessary
newTimeBin = 0.001; % in seconds
newTimeBase = (newTimeBin:newTimeBin:max(timestamps))';
Respiration_DS = interp1q(timestamps,Respiration,newTimeBase);
clear Respiration

% filter the data
order = 2; 
width = 21;
Respiration_DS_Filt = sgolayfilt(Respiration_DS,order,width);



%% load raw channels data

data_raw=[];
channelid = [1:channelNum];
for ch=1:channelNum 
    fname = sprintf('100_CH%d.continuous',channelid(ch));
    fprintf('loading ch %d [%d]\n',ch,channelid(ch));
    [data, timestamps, info] = load_open_ephys_data_faster(fullfile(myKsDir,fname));
   % data=data(1:20000); % cut down data
    data_raw(:,ch) = data;
end
disp('done');

%size(data_raw);

data_raw = data_raw.*info.header.bitVolts;
fs = info.header.sampleRate;

%% filter raw data 

% duration_to_analyze = 5*60*fs;
duration_to_analyze = size(data_raw, 1);

[b,a] = butter(1,[300 6000]/(fs/2),'bandpass'); 
data_bp = filtfilt(b,a,data_raw(1:duration_to_analyze,:));



%% make graphs
plotlim = 20000;
%dir = fullfile(myKsDir, strcat("cluster", num2str(cluster.id)));
%if exist(dir, 'dir')
%   rmdir(dir, 's')
%end
%mkdir(dir);
%figure();
sublplot(2,1,1)
plot(data_bp(1:plotlim,:)+100*repmat([1:32]',1,plotlim)','k');
%saveas(f, fullfile(dir, 'channels_and_ADC.png'));
sublplot(2,1,2)
plot(1:length(Respiration_DS),Respiration_DS_Filt);
    

