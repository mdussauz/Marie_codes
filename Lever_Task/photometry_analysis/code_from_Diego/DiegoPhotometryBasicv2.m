clear
clc
%addpath('C:\Users\deher\Dropbox\Smellocator project\Open Ephys\DH4\2020-02-12_19-13-54', 'C:\Users\deher\Dropbox\Functions');
%addpath('C:\Users\deher\Dropbox\Smellocator project\Open Ephys\DH1\2020-02-12_18-14-44', 'C:\Users\deher\Dropbox\Functions');

%% read and process event timestamps from open ephys file
[Events] = ParseOpenEphysEvents('all_channels.events');
% Channel 0: Trial start?
% Channel 1: Odor 1?
% Channel 2: Odor 2?
% Channel 3: Odor 3?
% Channel 5: water?
% StartTimeStamps      = [Events.Channel0.On Events.Channel0.Off];
% Odor1TimeStamps      = [Events.Channel1.On Events.Channel1.Off];
% Odor2TimeStamps      = [Events.Channel2.On Events.Channel2.Off];
% Odor3TimeStamps      = [Events.Channel3.On Events.Channel3.Off];
TrialTimeStamps      = [Events.Channel4.On Events.Channel4.Off];
% WaterTimeStamps      = [Events.Channel5.On Events.Channel5.Off];
% % delete first entry - just an empty trigger
% StartTimeStamps(1,:) = [];
% Odor1TimeStamps(1,:) = [];
% Odor2TimeStamps(1,:) = [];
% Odor3TimeStamps(1,:) = [];
TrialTimeStamps(1,:) = [];
% WaterTimeStamps(1,:) = [];

%% Read Photoreceivers and LEDs data from open ephys file

[modData_1, Timestamps_PR_1, info_PR_1]  = load_open_ephys_data('100_ADC3.continuous');
[modLED_1, Timestamps_LED_1, info_LED_1] = load_open_ephys_data('100_ADC5.continuous');

%% Constants

modFreq_1          = 211;       
modAmp_1           = 0.6;
samplingRate       = info_PR_1.header.sampleRate;
lowCutoff          = 15;

%% Prepare reference data and generate 90deg shifted reference data
   
shift_modLED_1     = modLED_1 - mean(modLED_1);                            % Remove DC offset
samplesPerPeriod   = (samplingRate/modFreq_1);
quarterPeriod      = round(samplesPerPeriod/4);
shift_modLED90_1   = circshift(shift_modLED_1,[1 quarterPeriod]);

%% Quadrature decoding and filtering                                       % Element-by-element array multiplication 
   
processedData0_1    = modData_1 .* shift_modLED_1;                         % 0 degrees data correction                          
processedData90_1  = modData_1 .* shift_modLED90_1;                        % 90 degrees data correction 

%% Low pass filter
    
norm_lowCutoff     = lowCutoff/(samplingRate/2);                           % CutOff normalized by half sampling rate 
[b, a]             = butter(5, norm_lowCutoff,'low');                      % '5th order' butterworth low pass filter

paddedData0_1        = processedData0_1(1:samplingRate,1);             
paddedData90_1       = processedData90_1(1:samplingRate,1);
demodDataFilt0_1     = filtfilt(b,a,[paddedData0_1; processedData0_1]);    % pad the data to suppress windows effect upon filtering
demodDataFilt90_1    = filtfilt(b,a,[paddedData90_1; processedData90_1]);        
processedData0f_1    = demodDataFilt0_1(samplingRate + 1: end, 1);
processedData90f_1   = demodDataFilt90_1(samplingRate + 1: end, 1);
 
demodData_1          = (processedData0f_1 .^2 + processedData90f_1 .^2) .^(1/2);

%% Correct for amplitude of reference

demodDataC_1         = demodData_1*(2/modAmp_1);
     
meanF0               = mean(demodDataC_1);                                 % Mean value of baseline
medianF0             = median(demodDataC_1);                               % Median value of baseline
pc1F0                = prctile(demodDataC_1,1);                            % 1% percentile value of baseline 
pc5F0                = prctile(demodDataC_1,5);                            % 5% percentile value of baseline
pc10F0               = prctile(demodDataC_1,10);                           % 10% percentile value of baseline 
pc20F0               = prctile(demodDataC_1,20);                           % 20% percentile value of baseline
pc40F0               = prctile(demodDataC_1,40);                           % 40% percentile value of baseline 
pc80F0               = prctile(demodDataC_1,80);                           % 80% percentile value of baseline


%% Read the stimulus list text file

StimList = repmat((1:20),1,5);


    %% for each stimulus type

    for i = 1:numel(unique(StimList))
  
        reps            = find(StimList==StimList(i));                     % get indices of all repeats
        DFF_1           = [];                                              % creates array for DFF_1
        
        for j = 1:numel(reps)                                              % for each trial
        
        %% trial info
        
        idxON    = find(Timestamps_PR_1==TrialTimeStamps(reps(j),1),1,'first');
        idxOFF   = find(Timestamps_PR_1==TrialTimeStamps(reps(j),2),1,'first');
        nidx     = length(idxON:idxOFF);
        myTrace  = demodDataC_1((idxON-2500):idxOFF);    
        
        %% F0 option
        
        F0       = pc1F0;
        
        %% dF/F0
        
        DFF                    = 100*((myTrace-F0)/F0);
        DFF_1(j,1:length(DFF)) = DFF;
        DFF_1(DFF_1==0)        = NaN;                                      % replace 0 values by NaN's
        
        end
      
    subplot(4,5,i);
    plot(DFF_1');
    
    end

% rescale plots
for i = 1:10
    subplot(4,5,i);
    set(gca,'YLim',[-5 30]);
    set(gca,'XLim',[0 200000]);
    trial_a = vline(25000,'r');
end
for i = 11:15
    subplot(4,5,i);
    set(gca,'YLim',[-5 30]);
    set(gca,'XLim',[0 200000]);
    trial_b = vline(25000,'r');
end
for i = 16:20
    subplot(4,5,i);
    set(gca,'YLim',[-5 30]);
    set(gca,'XLim',[0 200000]);
    trial_c = vline(25000,'r');
end

% figure(13)
% plot(modLED_1)
% set(gca, 'XLim',[0 100000000]);                                               % Limit x axis range
% set(gca, 'YLim',[0 0.7]);
% title('Whole Session DH4')
% xlabel('sample')
% ylabel('modulated LED (V)')

%F0 baseline values histogram

% [H1,P1,STATS1]  = chi2gof(demodDataC_1)                                    % normal fit statistics
% [H2,P2,STATS2]  = chi2gof(log(demodDataC_1))                               % log-normal fit statistics
% 
figure(2);
fdistnormal = fitdist(demodDataC_1,'normal');
histfit(demodDataC_1,100,'normal');   
% set(gca, 'XLim',[0.3 1.7]);                                               % Limit x axis range
% set(gca, 'YLim',[0 14000000]);
% xlabel('demodulated signal (a.i.)')
% ylabel('frequency')
% title('Normal distribution')
% mean_norm = vline(meanF0,'b');
% median_norm = vline(medianF0,'g');
% fit_mean_norm = vline(fdistnormal.mu,'p');
% pc1 = vline(pc1F0,'r');
% text(0.4, 2000000, 'percentile 1 = 0.8681');
% text(0.4, 2500000, 'mean = 0.9631');
% text(0.4, 3000000, 'median = 0.9660');

%text(0.4, 50000, 'Normal chi^2 p < 0.05')  
% 
% figure(3);
% fdistlognormal  = fitdist(demodDataC_1,'lognormal');
% histfit(demodDataC_1,100,'lognormal');   
% xlabel('fluorescence')
% ylabel('frequency')
% title('Log-normal distribution')
% mean_log = vline(meanF0,'b');
% median_log = vline(medianF0,'g');
% fit_mean_lognorm = vline(exp((fdistlognormal.mu)+((fdistlognormal.sigma^2)/2)),'p');
% %text(0.4, 500, 'Log-normal chi^2 p < 0.05');
% 
% figure(4);
% histogram(demodDataC_1)
% xlabel('Demodulated signal (a.i.)')
% ylabel('frequency')
% pc1 = vline(pc1F0,'r');
% title('Percentile 1');
% text(0.4, 300000, 'percentile 1 = 0.8681');
% 
% figure(5);
% histogram(demodDataC_1)
% xlabel('fluorescence')
% ylabel('frequency')
% pc40 = vline(pc20F0,'r');
% title('Percentile 20');
% text(0.4, 300000, 'percentile 40 = 0.9267');
% 
% figure(6);
% histogram(demodDataC_1)
% xlabel('fluorescence')
% ylabel('frequency')
% pc40 = vline(pc40F0,'r');
% title('Percentile 40');
% text(0.4, 300000, 'percentile 40 = 0.9552');
% 
% figure(7);
% histogram(demodDataC_1)
% xlabel('fluorescence')
% ylabel('frequency')
% pc40 = vline(pc80F0,'r');
% title('Percentile 80');
% text(0.4, 300000, 'percentile 80 = 0.9992');