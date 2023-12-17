%AnalyzingDriftwithDriftMap

%% paths

% Free Lever
SessionName = 'S12/2023-08-04_14-29-23';

if strcmp(computer, 'MACI64')
    ephyspath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Ephys/';
else
    ephyspath = '/mnt/data/Processed/Ephys/';
end

myKsDir = fullfile(ephyspath,SessionName);

%% Get relevant information: 
[spikeTimes, spikeAmps, spikeDepths, spikeSites] = ksDriftmap(myKsDir);

%% Plot drift map: 
figure; plotDriftmap(spikeTimes, spikeAmps, spikeDepths); 
% plot spike depth by spiketimes 
% the depth of a spike is defined as the center of mass of the features
% (1st PC component)
% is sum(coords.*features)/sum(features)
% every 100 is basically a different channel 

%% Computing useful properties of spikes and templates
% sp = loadKSdir(myKsDir);
% 
% [spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
%     templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
% 
% [spikeWidths, tempWidths, clusterWidths] = computeSpikeWidths(tempsUnW, sp.spikeTemplates);

