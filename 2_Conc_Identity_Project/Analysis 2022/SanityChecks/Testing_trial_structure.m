% Testing_trial_structure 

SessionInfo = ReadSessionDatatable(); 

if strcmp(computer,'PCWIN64')
    KSpath = 'Z:\mdussauz\ephysdata\Conc_id_exp';
else
    KSpath = '/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp';
    addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
end

i = 16; % session number to check 
mousename = SessionInfo.MouseName(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
stimfilename = SessionInfo.StimFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
filenameKS = fullfile(KSpath,mousename{i}, myKsDir{i},'Record Node 106');

%bypass = '/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp/E2/2022-06-11_17-29-14/Record Node 106/all_channels.events';
%bypass = '/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp/E2/2022-06-11_13-57-38/Record Node 106/all_channels.events';
%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Testing strategies to solve potential extra TTLs problem 
strategy = 1; % 2 is the old stratgey 
switch strategy
    case 1 %trying new strategy 
        [StimTime, StimList, Nrepeats] = ReadStimFile([stimfilename{i}, '.txt']);

    case 2 % this strategy works if there are no extra errors before the actual exp begins
            %% Get Events and correct for ephys offset
            [Events] = ParseOpenEphysEvents(fullfile(filenameKS,'all_channels.events'));
            %[Events] = ParseOpenEphysEvents(bypass);
            % trial events are in Channel0, odor events are in Channel1

            TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
            TrialTimestamps(2,:) = []; % delete imaging empty trigger
            OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

            %     TrialTimestamps = TrialTimestamps - offset; %adjust for offset
            %     OdorTimestamps = OdorTimestamps - offset; % adjust for offset

            %% Conclusion
            % Odortimestamps(1) happens before trial timestamps(1,1)
            % All E2 and E3 sessions have on extra odor on/off and trial on/off when
            % the DAC resets (less than a sec)
            % as well as an extra trial on/off for imaging with length = start pre stim
            % to end post stim
            if OdorTimestamps(1,2) - OdorTimestamps(1,1) < 2
                TrialTimestamps(1,:) = [];
                OdorTimestamps(1,:) = [];
            end

            Trial_length = TrialTimestamps(:,2) - TrialTimestamps(:,1);
            Odor_length = OdorTimestamps(:,2) - OdorTimestamps(:,1);

end