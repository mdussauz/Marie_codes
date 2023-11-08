function  [StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename)

%%
%written by MD
%read the stimulus txt file and get the settings:
% -StimTime(1) = pre stimulus time
% -StimTime(2) = stimulus time
% -StimTime(3) = post stimulus time
% -StimTime(4) = ITI
% -StimList = list and order of odors during session
% - Nrepeats = number of repeats for each odor 

%%

if strcmp(computer, 'PCWIN64')
    addpath(genpath('Z:\mdussauz\PhotonCerber_Stimuli_on_server'))
else
    addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
end

%stimfilename = '190910_17_01.txt';

fileID = fopen(stimfilename,'r');
Stimulus = fscanf(fileID,'%f');
StimTime = Stimulus(2:5);
StimList = Stimulus(7:end);
Nrepeats = Stimulus(6);

end

