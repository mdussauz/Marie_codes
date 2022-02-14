function  [StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename)

%%
%written by MD
% Inputs = 
%stimfilename = 'path/190910_17_01.txt';
%read the stimulus txt file and get the settings:
% Outputs =
% -StimTime(1) = pre stimulus time
% -StimTime(2) = stimulus time
% -StimTime(3) = post stimulus time
% -StimTime(4) = ITI
% -StimList = list and order of odors during session
% - Nrepeats = number of repeats for each odor 

%%

%addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
%stimfilename = '190910_17_01.txt';

fileID = fopen(stimfilename,'r');
Stimulus = fscanf(fileID,'%f');
StimTime = Stimulus(2:5);
StimList = Stimulus(7:end);
Nrepeats = Stimulus(6);

end

