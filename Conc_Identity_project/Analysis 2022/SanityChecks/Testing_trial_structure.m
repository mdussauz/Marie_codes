% Testing_trial_structure 

SessionInfo = ReadSessionDatatable(); 
KSpath = 'Z:\mdussauz\ephysdata\Conc_id_exp'; 
myKsDir = SessionInfo.RecFile(strcmp( 'Y',SessionInfo.IncludeInAnalysis) & strcmp('AON',SessionInfo.BrainRegion));
filenameKS = fullfile(KSpath, myKsDir{2},'Record Node 106');

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 
%% Offset
[offset] = AdjustClockOffset(filenameKS);
offset = offset/sampleRate;

%% Get Events and correct for ephys offset
    [Events] = ParseOpenEphysEvents(fullfile(filenameKS,'all_channels.events'));
    % trial events are in Channel0, odor events are in Channel1
    
    TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
    TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
    OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
    

    
    TrialTimestamps = TrialTimestamps - offset; %adjust for offset
    OdorTimestamps = OdorTimestamps - offset; % adjust for offset
    
%% Conclusion
% Odortimestamps(1) happens before trial timestamps(1,1)
if OdorTimestamps(1,2) - OdorTimestamps(1,1) < 2
    TrialTimestamps(1,:) = [];
    OdorTimestamps(1,:) = [];
end