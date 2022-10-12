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

% make the spike structure 
[allcluster] = MakeSpikeStucture('AON', 'all');

%% If exist, load spike structure 