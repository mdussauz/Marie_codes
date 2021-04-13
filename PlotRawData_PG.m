%plot raw data 

%% OVERVIEW
%{

1 - load some data from one tetrode 
2 - remove very low  and very high frequency data, leaving only the freq.
band where spikes live
3 - apply a threshold to find the times of these spikes
4- for each of these spikes, remember the max. amplitude, and the waveform
5 - look at the clustering of amplitudes across tetrode channels

%}


%% 1- load raw data

set(0,'DefaultFigureWindowStyle','docked'); % fix matlab's figure positioning bug

% raw data available on
% https://drive.google.com/drive/folders/1CwFcErgp3F3D6I2TB_hTtW1JAQB21TAC?usp=sharing
%
%datapath='C:\Users\mdussauz\Desktop\Spike_Sorting\TENSS\spike_sorting_tutorial_Nacho\tenss_example_data'; % <- edit this
%datapath='/mnt/data/N5/2019-09-24_15-13-13'
datapath='/Users/mdussauz/Desktop/Analysis/M04/2021-03-30_16-46-03';

data_raw=[];
channelid = [17:32];
nb_tetrodes = numel(channelid)/4;

colors = linspecer(nb_tetrodes ,'qualitative');

i =1;

%colors =  ['g';'b';'pd';'pl';'o']
%Plot_Colors(colors(1))

figure; 
for ch=1:numel (channelid) 
    fname = sprintf('100_CH%d.continuous',channelid(ch));
    fprintf('loading ch %d/32 [%d]\n',ch,channelid(ch));
    [data_raw, timestamps, info] = load_open_ephys_data_faster(fullfile(datapath,fname));
   % data=data(1:20000); % cut down data
    %data_raw(:,ch) = data;
    %data_raw = data;
    %size(data_raw);

    data_raw = data_raw.*info.header.bitVolts;
    fs = info.header.sampleRate;
    
    plotlim = [10:20];
    time_to_plot = [(plotlim(1)*fs):(plotlim(end)*fs)];

    % duration_to_anlyse = 5*60*fs;
    duration_to_anlyse = size(data_raw, 1);

    [b,a] = butter(1,[300 6000]/(fs/2),'bandpass'); 
    data_bp = filtfilt(b,a,data_raw(1:duration_to_anlyse,:));
    
    clear data_raw
    
    if ismember(ch,[1;2,;3;4])
        i =1;
    elseif ismember(ch,[5;6,;7;8])
        i =2;
    elseif ismember(ch,[9;10;11;12])
        i = 3;
    elseif ismember(ch,[13;14;15;16])
        i = 4;
    end
        
    plot(data_bp(time_to_plot,:)+100*ch,'Color', colors(i,:)); hold on

end
% disp('done');
% 
% %%
% size(data_raw);
% 
% data_raw = data_raw.*info.header.bitVolts;
% fs = info.header.sampleRate;

%% filter data

% plotlim=10000;
% 
% % duration_to_anlyse = 5*60*fs;
% duration_to_anlyse = size(data_raw, 1);
% 
% [b,a] = butter(1,[300 6000]/(fs/2),'bandpass'); 
% data_bp = filtfilt(b,a,data_raw(1:duration_to_anlyse,:));
% 
% figure(3); clf;
% plot(data_bp(1:plotlim,:)+100*repmat([1:numel(channelid)]',1,plotlim)','k');

