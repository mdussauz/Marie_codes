% MasterFile
% -- written by MD
% comment sections based on needs 

% WIP - note to self 
% include part for analysing 
% taking into account old trial strc for previous round of exp 
% taking into account structure of Diego's photometry data

%% -- COMPILE the spike structure 
% comment this section if structure already exists and instead run next
% section

% extract info about session
%SessionInfo = ReadSessionDatatable(FilePath, FileName); 
SessionInfo = ReadSessionDatatable(); 

% make the spike structure - specify AON or APC
% allclusters contains both conc and id exp 

%[allclusters] = MakeSpikeStucture(SessionInfo, 'AON', 'all');
[allclusters] = MakeSpikeStucture(SessionInfo, 'APC', 'all');

%% If exist, load spike structure 

%% Compute PSTH 5D for dPCA
[goodcluster,spikes] = ComputePSTHMultiD_new(allclusters, "Conc",5); 
%[goodcluster,PSTH] = ComputePSTHMultiD_new(allclusters, "Id",4); 

%% Smooth PSTH
[PSTH] = SmoothPSTH_new(spikes, 100, "Conc");
%[smoothPSTH] = SmoothPSTH_new(PSTH, 100, "Id");

%% Plotting Conc responses
 plotConcResponseMatrix_new(smoothPSTH);

%% Run dPCA
[W] = RundPCA(smoothPSTH);