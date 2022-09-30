function [allcluster] = MakeSpikeStucture(BrainRegion, Mouse)
% -- written by MD
% -- create structure array for all neurons all sessions all mice
% user must specify what to extract: 
% sessions coming from AON or APC mouse 
% whether to use all mice or just a specific mouse 
% whether to build it for a specific session   

%% defaults
sampleRate = 30000; % Open Ephys acquisition rate = 30 kHz 

%% Stimulus File

SessionInfo.RecFile(strcmp( 'N',SessionInfo.IncludeInAnalysis) & strcmp( 'E2',SessionInfo.MouseName))