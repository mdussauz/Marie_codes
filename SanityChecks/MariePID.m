Exp = ''; % leave empty (i.e. '') if you want to select file in folder

%% Define datapath based on the computer being used
if strcmp(computer, 'PCWIN64') % Marie's remote desktop
    KSpath = 'Z:\mdussauz\PID_test';
    StimFolder = 'Z:\mdussauz\PhotonCerber_Stimuli_on_server';
    
else % Marie's work linux machine
    KSpath = '/mnt/grid-hs/mdussauz/PID_test';
    StimFolder = 'mnt/grid-hs/mdussauz/PhotonCerber_Stimuli_on_server';
end

%%
switch Exp
    case ''
        [~,recsessionpath] = uigetfile([KSpath, '\.continuous'], 'Select PID signal file');
        [stimefilename, StimFolder] = uigetfile([StimFolder,'\.txt'], 'Select stimulus file');
        stimfilepath = fullfile(StimFolder,stimefilename);
        
    
    %PID experiments on stripped machine which hasn't encountered odor yet:
    %NB - for the Blank and oil runs, only 1 rep with stim non randomized
    case '220525_Conc_Blank'
        recsessionpath = fullfile(KSpath, '2022-05-25_08-58-10_Conc_Blank\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220525\220525_10_20.txt');

    case '220525_Id_Blank'
        recsessionpath = fullfile(KSpath, '2022-05-25_08-58-10_16Odor_Blank\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220525\220525_11_45.txt');

    case '220525_Conc_Oil'
        recsessionpath = fullfile(KSpath, '2022-05-25_08-58-10_Conc_Oil\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220525\220525_13_51.txt');

    case '220525_Id_Oil'
        recsessionpath = fullfile(KSpath, '2022-05-25_08-58-10_16Odor_Oil\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220525\220525_13_28.txt');
    
    %PID experiments on stripped machine
    case '220530_Conc'
        recsessionpath = fullfile(KSpath, '2022-05-30_18-19-56\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220530\220530_18_20.txt');

    case '220531_Id'
        recsessionpath = fullfile(KSpath,'2022-05-31_08-54-29_16Odor\Record Node 104');
        stimfilepath = fullfile(StimFolder,'220531\220531_9_02.txt');
        
   %PID experiments after 2nd round of conc/id experiments -'dirty' machine
end
    
%% read the stimulus list text file
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilepath);

baseline = 1:StimTime(1); % pre stimulus time
TrialLength = StimTime(1)+StimTime(2)+StimTime(3);

StimList(find(StimList==0))= []; %zeros are labview errors
Nodors = numel(unique(StimList));% number of odors 
NTrials = numel(StimList);
%% read PID data from open ephys file
if  strcmp(Exp,'220525_Conc_Oil')% this session has a naming issue
    [PID, Timestamps, info] = load_open_ephys_data(fullfile(recsessionpath,'100_1_2.continuous'));
elseif strcmp(Exp,'220525_Id_Oil')
    [PID, Timestamps, info] = load_open_ephys_data(fullfile(recsessionpath,'100_1_2.continuous'));
else
    [PID, Timestamps, info] = load_open_ephys_data(fullfile(recsessionpath,'100_1.continuous'));
end 

% adjust timestamps to account for the start offset in OEPS
OepsSampleRate = 30000; % Open Ephys acquisition rate

%% read and process event timestamps from open ephys file
if  strcmp(Exp,'220525_Conc_Oil')% this session has a naming issue
    [Events] = ParseOpenEphysEvents(fullfile(recsessionpath,'all_channels_2.events'));
elseif strcmp(Exp,'220525_Id_Oil')
    [Events] = ParseOpenEphysEvents(fullfile(recsessionpath,'all_channels_2.events'));
else 
    [Events] = ParseOpenEphysEvents(fullfile(recsessionpath,'all_channels.events'));
end
% trial events are in Channel1, odor events are in Channel2

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
if length(TrialTimestamps)>NTrials
    TrialTimestamps = TrialTimestamps((end - NTrials+1):end,:);
end

%% Make grid of subplots
if Nodors == 20
    gridX = 4;
    gridY = 5;
elseif Nodors == 16
    gridX = 4;
    gridY = 4;
end
SamplePerSec = OepsSampleRate/1000;
%% for each stimulus type
for i = 1:Nodors
    % get indices of all repeats
    reps = find(StimList==StimList(i));
    traces = [];
    % for each repeat
    for j = 1:numel(reps)
        idxON = find(Timestamps==TrialTimestamps(reps(j),1),1,'first');
        idxOFF = find(Timestamps==TrialTimestamps(reps(j),2),1,'first');
        nidx = length(idxON:idxOFF);
        mytrace = interp1(1:nidx,PID(idxON:idxOFF),SamplePerSec:SamplePerSec:nidx);
        %mytrace = PID(idxON:idxOFF);
        % subtract baseline
        mytrace = mytrace - mean(mytrace(baseline));
        traces(j,1:length(mytrace)) = mytrace;
    end
    subplot(gridX,gridY,StimList(i));
    plot(traces');
end

