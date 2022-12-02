% MasterFile
% -- written by MD
% comment sections based on needs 

% WIP - note to self 
% taking into account old trial strc for previous round of exp 
% taking into account structure of Diego's photometry data

%% Specify what region and what mouse (all or specific one)
BrainRegion = 'AON'; % choices: AON or APC
Mouse = 'E3'; % choices: all, E2, E3, or E6
WhatExp = 'Conc'; %choices: Conc or Id

%% -- Extract info about recording sessions
%SessionInfo = ReadSessionDatatable(FilePath, FileName); 
SessionInfo = ReadSessionDatatable(); 

%% -- If exist, load spike structure 
if strcmp(computer, 'PCWIN64') % Marie's remote desktop 
    datapath = 'C:\Users\Marie\Documents\data\allclusters\2022';
else % Marie's work linux machine
    datapath = '/mnt/data/allclusters/2022';
end 

savefile = ['allclusters_',WhatExp, '_',BrainRegion,'_',Mouse,'.mat'];
savepath = fullfile(datapath,savefile);

if exist(savepath, 'file')
    reply = input('This spike structure already exists. \nDo you want to load it? Y/N [Y]: ','s');
    if strcmp(reply,'Y')
        load(savepath)
    end
end

%% -- COMPILE the spike structure if doesn't exist or user wants to overwrite
% make the spike structure - specify AON or APC and which mouse 
% allclusters contains both conc and id exp 
if ~exist(savepath, 'file') || strcmp(reply,'N') 
    [allclusters] = MakeSpikeStucture(SessionInfo, BrainRegion, Mouse);
    save(savepath, 'allclusters')
end 

%% Compute PSTH 5D for dPCA
% 1st output is str array for the chosen exp type
% 2nd output is the spikes (0 or 1) matrix

% 3rd input is optional and allows to specify whether to make a 4 or 5 dim
% matrix
[goodcluster,spikes] = ComputePSTHMultiD_v2(allclusters, WhatExp); 

%% Smooth PSTH
[smoothPSTH] = SmoothPSTH_v2(spikes, 100, WhatExp);

%% Plotting Conc responses
plotConcResponseMatrix_v2(smoothPSTH);

%% Run dPCA
%[W, explVar] = RundPCA(smoothPSTH);
[W, explVar] = RundPCA(smoothPSTH(:,1:4,:,:,:));

%% Plot dPCA plots
PlotdPCA(smoothPSTH(:,1:4,:,:,:), W, explVar)

%% Run odor generalizer decoder
OdorGeneralizer(smoothPSTH)

%% Run large odor set decoder 
LargeOdorSetDecoder(smoothPSTH)