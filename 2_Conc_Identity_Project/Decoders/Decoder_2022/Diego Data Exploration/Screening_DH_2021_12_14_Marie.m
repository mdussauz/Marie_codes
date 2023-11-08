close all
clearvars

%% Add path

addpath(genpath('C:\Users\deher\Dropbox\APC feedback\Functions')); % Functions: vline.m 

%% Load data folder        

current_dir     = 'D:\Cancer_project\Cancer_project\APC\Post_training'; % Variable con el path a la carpeta con las imagenes

cd(current_dir)                     % define 'mouse training sessions directory'

dataset         = dir;              % variable with list of the sessions in the folder
clear_index     = find([dataset.isdir]);

dataset(clear_index,:)  = [];

file_num        = size(dataset,1);  % Number of files/elements in dataset folder

Data_comp       = NaN(file_num,92720,16,10);

for nFile = 1:file_num   
    
    session_name    = dataset(nFile).name;  % define dataset to load
    load(session_name)
    
    %% Session properties
    
    num_Trials          = SessionData.nTrials; 
    type_Trials         = SessionData.TrialTypes';
    duration            = SessionData.TrialSettings(1).GUI.NidaqDuration;
    sampleRate          = SessionData.TrialSettings(1).GUI.NidaqSamplingRate;
    baseline_begin      = SessionData.TrialSettings(1).GUI.BaselineBegin;
    baseline_end        = SessionData.TrialSettings(1).GUI.BaselineEnd;
    modFreq             = SessionData.TrialSettings(1).GUI.LED.LED1(3);         % Modulation frequency
    modAmp              = SessionData.TrialSettings(1).GUI.LED.LED1(2);         % Modulation amplitude
    
    decimateFactor      = 1;
    TimeToZero          = 10;
    
    %% Data demodulation
    
    for i = 1:num_Trials                                                          % Trial 'for-loop' / i: trial number for the loop

        Photoreceiver_Data  = SessionData.Photoreceiver1Data{1,i};                  % 1: Modulated LED waveform 2: Modulated data acquired   
        Length_PD           = length(Photoreceiver_Data);                           % Number of samples Photoreceiver_Data matrix
        refData             = Photoreceiver_Data(:,2);                              % Modulated LED waveform
        rawData             = Photoreceiver_Data(:,1);                              % Acquired modulated data
        
        refData             = refData(1:length(rawData),1);   % match length of refData and rawData
        refData             = refData - mean(refData);          % remove DC offset
        samplesPerPeriod    = (1/modFreq)/(1/sampleRate);
        quarterPeriod       = round(samplesPerPeriod/4);
        refData90           = circshift(refData,[1 quarterPeriod]);

        %% Quadrature decoding and filtering
        processedData_0     = rawData .* refData;
        processedData_90    = rawData .* refData90;

        %% Filter
        % For filtering
        lowCutoff           = 15;
        pad                 = 1;
    
        
        lowCutoff = lowCutoff/(sampleRate/2); % normalized CutOff by half SampRate (see doc)
        [b, a] = butter(5, lowCutoff, 'low'); 

        % pad the data to suppress windows effect upon filtering
        if pad == 1
            paddedData_0        = processedData_0(1:sampleRate, 1);
            paddedData_90       = processedData_0(1:sampleRate, 1);
            demodDataFilt_0     = filtfilt(b,a,[paddedData_0; processedData_0]);
            demodDataFilt_90    = filtfilt(b,a,[paddedData_90; processedData_90]);        
            processedData_0     = demodDataFilt_0(sampleRate + 1: end, 1);
            processedData_90    = demodDataFilt_90(sampleRate + 1: end, 1);
        else
            processedData_0     = filtfilt(b,a,processedData_0);
            processedData_90    = filtfilt(b,a,processedData_90); 
        end

        demodData = (processedData_0 .^2 + processedData_90 .^2) .^(1/2);

        %% Correct for amplitude of reference
        demodData = demodData * 2/modAmp;

        %% Expected Data set
        SampRate                    = sampleRate/decimateFactor;
        ExpectedSize                = duration*SampRate;
        Data                        = NaN(ExpectedSize,1);
        TempData                    = decimate(demodData,decimateFactor);
        Data(1:length(TempData))    = TempData;

        %% DF/F calculation
        Fbaseline = mean(Data(baseline_begin*SampRate:baseline_end*SampRate));
        DFF       = 100*(Data-Fbaseline)/Fbaseline;

        %% Time
        Time        = linspace(0,duration,ExpectedSize);
        Time        = Time' - TimeToZero;

        %% Raw Data
        ExpectedSizeRaw             = duration*sampleRate;
        DataRaw                     = NaN(ExpectedSizeRaw,1);
        DataRaw(1:length(rawData))  = rawData;

        TimeRaw = linspace(0,duration,ExpectedSizeRaw);
        TimeRaw = TimeRaw' - TimeToZero;

        %% NewDataSet
        Data_Demod(:,1)     = Time;
        Data_Demod(:,2)     = Data;
        Data_Demod(:,3)     = DFF;

        Data_Raw(:,2)       = DataRaw;
       
        DFF_session(:,i)    = DFF;
    end
    
    clearvars -except current_dir DFF_session Time nFile file_num ...
                      num_Trials SessionData type_Trials dataset ...
                      decimateFactor i Data_comp
    
    %% Trial type index and data parsing
    
    for u = 1:16
       
        odor_index              = find(type_Trials == u);
        
        Data_comp(nFile,:,u,1:size(odor_index,1))  = DFF_session(:,odor_index);

    end
    
    clearvars u
end

    time_a = 1; % -10
    time_b = 30500/decimateFactor;  % -5
    time_c = 61000/decimateFactor; % 0
    time_d = 91500/decimateFactor; % 5  
    time_e = 63135/decimateFactor;  % 0.35
    
    odor_list = {'isoamylamine','ehtyl tiglate','heptanal','allyl tiglate', ...
                 'ethyl butyrate','ethyl propionate','isoamyl acetate',     ...
                 'ethyl valerate','1-butanol','allyl butyrate',             ...
                 'propyl butyrate','ethyl heptanoate','hexanal','limonene', ...
                 'oil','air'};   

    Data_compilation.Data           = Data_comp;
    Data_compilation.odor_list      = odor_list';
    Data_compilation.events(1).Legend    = {'-10 s','-5 s','0 s','0.35 s','5 s'};
    Data_compilation.events(1).Time      = [-10 -5 0 0.35 5];
    Data_compilation.events(1).Frame     = [time_a time_b time_c time_e time_d];
   
