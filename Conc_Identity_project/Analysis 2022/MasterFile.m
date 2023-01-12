% MasterFile
% -- written by MD
% comment sections based on needs 

% WIP - note to self 
% taking into account old trial strc for previous round of exp 
% taking into account structure of Diego's photometry data

%% Specify what region and what mouse (all or specific one)
BrainRegion = 'APC'; % choices: AON or APC
Mouse = 'E6'; % choices: all, E2, E3, or E6
WhatExp = 'Conc'; %choices: Conc or Id

%% -- Extract info about recording sessions
%SessionInfo = ReadSessionDatatable(FilePath, FileName); 
SessionInfo = ReadSessionDatatable(); 

%% -- If exist, load spike structure 
% -- COMPILE the spike structure if doesn't exist or user wants to overwrite
if strcmp(computer, 'PCWIN64') % Marie's remote desktop 
    datapath = 'C:\Users\Marie\Documents\data\allclusters\2022';
else % Marie's work linux machine
    datapath = '/mnt/data/allclusters/2022';
end 

savefile = ['allclusters_',BrainRegion,'_',Mouse,'.mat'];
savepath = fullfile(datapath,savefile);

if exist(savepath, 'file')
    reply = input('This spike structure already exists. \nDo you want to load it? Y/N [Y]: ','s');
    if strcmp(reply,'Y')
        load(savepath)
    end
end

% make the spike structure - specify AON or APC and which mouse 
% allclusters contains both conc and id exp 
if ~exist(savepath, 'file') || strcmp(reply,'N') 
    [allclusters] = MakeSpikeStucture(SessionInfo, BrainRegion, Mouse);
    save(savepath, 'allclusters') %save structure
end 

%% Compute PSTH 5D for dPCA
% 1st output is str array for the chosen exp type
% 2nd output is the spikes (0 or 1) matrix

% 3rd input is optional and allows to specify whether to make a 4 or 5 dim
% matrix
[goodcluster,spikes] = ComputePSTHMultiD_v2(allclusters, WhatExp); 

%% Smooth PSTH from the spikes matrix - in Hz 
[smoothPSTH] = SmoothPSTH_v2(spikes, 100, WhatExp);

%% Plotting Conc responses
plotConcResponseMatrix_v2(smoothPSTH);

%% Plotting Conc z scored responses
plotZScoredOdorConcMatrix_v3(smoothPSTH);

%% Plotting Conc responses with respect to oil
plotOilZScoredOdorConcMatrix(smoothPSTH)
%% Run dPCA
%[W, explVar, whichMarg] = RundPCA(smoothPSTH); %with oil and all conc
%[W, explVar, whichMarg] = RundPCA(smoothPSTH(:,1:4,:,:,:)); %without oil and all conc
[W, explVar, whichMarg] = RundPCA(smoothPSTH(:,1:4,2:4,:,:)); %without oil and 3 highest conc

%% Plot dPCA plots
%PlotdPCA(smoothPSTH, W, explVar, whichMarg) %with oil and all conc
%PlotdPCA(smoothPSTH(:,1:4,:,:,:), W, explVar, whichMarg) %without oil and all conc
PlotdPCA(smoothPSTH(:,1:4,2:4,:,:), W, explVar, whichMarg) %without oil and 3 highest conc

%% Run odor generalizer decoder - generalization to novel concentration
[GENERAL_PERF] = OdorGeneralizer(smoothPSTH);

%% Run large odor set decoder 
LargeOdorSetDecoder(smoothPSTH)

%% Run Concentration Invariant Classifier(smoothPSTH)
[incasefml] = ConcInvariantClassifier(smoothPSTH);

%% Run identity and concentrtaion calling classifier 
[GENERAL_PERF] = IdentityAndConcentrationCalling(smoothPSTH);

%% Run Concentration predictor 
[GENERAL_PERF_Pre] = ConcentrationPredictor(smoothPSTH);

