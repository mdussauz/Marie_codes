%ReplayStretchResponses
%Script

%% paths
%SessionName = fullfile('S12','S12_20230727_r0_processed.mat');
SessionName = fullfile('O3','O3_20211005_r0_processed.mat');
ChosenUnits = [];

if strcmp(computer, 'MACI64')
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    ProcessedBehaviorPath = '/mnt/data/Processed/Behavior/';
end

MySession = fullfile(ProcessedBehaviorPath,SessionName);
%% Plot the selected units
PlotReplayResponsesMD(SessionName,ChosenUnits);