Exp = 'Conc';

%% Define datapath depending on the computer being used
if strcmp(computer, 'PCWIN64') % Marie's remote desktop 
    addpath(genpath('Z:\mdussauz\PhotonCerber_Stimuli_on_server'))
    addpath(genpath('C:\Users\Marie\Documents\data\PhotonCerber_Stimuli_on_server'))
    KSpath = 'Z:\mdussauz\PID_test'; 
    
else % Marie's work linux machine
    addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
    addpath(genpath('mnt/grid-hs/mdussauz/PhotonCerber_Stimuli_on_server'))
    KSpath = '/mnt/grid-hs/mdussauz/PID_test';

end

%%
switch Exp
    case '220530_Conc'
        addpath(genpath('C:\Users\Marie\Documents\data\Final PID test - May 2022\2022-05-30_18-19-56\Record Node 104'))
        stimfilename = 'C:\Users\Marie\Documents\data\Final PID test - May 2022\2022-05-30_18-19-56\Record Node 104\220530_18_20.txt';
        gridX = 4;
        gridY = 5;
    case '220531_Id'
        addpath(genpath('C:\Users\Marie\Documents\data\Final PID test - May 2022\2022-05-31_08-54-29_16Odor\Record Node 104'))
        stimfilename = 'C:\Users\Marie\Documents\data\Final PID test - May 2022\2022-05-31_08-54-29_16Odor\Record Node 104\220531_9_02.txt';
        gridX = 4;
        gridY = 4;

    case '220525_Conc_Blank'
        2022-05-25_08-58-10_Conc_Blank

    case '220525_Id_Blank'

        2022-05-25_08-58-10_16Odor_Blank

    case '220525_Conc_Oil'

        2022-05-25_08-58-10_Conc_Oil

    case '220525_Id_Oil'
        2022-05-25_08-58-10_16Odor_Oil
    
end

%% read PID data from open ephys file
[PID, Timestamps, info] = load_open_ephys_data('100_1.continuous');

% adjust timestamps to account for the start offset in OEPS
OepsSampleRate = 30000; % Open Ephys acquisition rate

%% read and process event timestamps from open ephys file
[Events] = ParseOpenEphysEvents('all_channels.events');
% trial events are in Channel1, odor events are in Channel2

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
% delete first entry - just an empty trigger
TrialTimestamps(1:2,:) = [];
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
OdorTimestamps(1,:) = [];

%% read the stimulus list text file
%StimList = repmat((1:20),1,5);
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
baseline = 1:4000;
%baseline = 1:StimTime(1); % pre stimulus time

%% for each stimulus type
for i = 1:numel(unique(StimList))
    % get indices of all repeats
    reps = find(StimList==StimList(i));
    traces = [];
    % for each repeat
    for j = 1:numel(reps)
        idxON = find(Timestamps==TrialTimestamps(reps(j),1),1,'first');
        idxOFF = find(Timestamps==TrialTimestamps(reps(j),2),1,'first');
        nidx = length(idxON:idxOFF);
        %mytrace = interp1(1:nidx,PID(idxON:idxOFF),32:32:nidx);
        mytrace = PID(idxON:idxOFF);
        % subtract baseline
        %mytrace = mytrace - mean(mytrace(baseline));
        traces(j,1:length(mytrace)) = mytrace;
    end
    subplot(gridX,gridY,StimList(i));
    plot(traces');
end

% rescale plots
% for i = 1:10
%     subplot(4,5,i);
%     set(gca,'YLim',[-0.2 0.8]);    
% end
% for i = 11:15
%     subplot(4,5,i);
%     set(gca,'YLim',[-0.2 0.8]);    
% end
% for i = 16:20
%     subplot(4,5,i);
%     set(gca,'YLim',[-0.2 0.4]);    
% end