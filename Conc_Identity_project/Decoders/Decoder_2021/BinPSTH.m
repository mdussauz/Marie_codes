
function [binnedPSTH] = BinPSTH(PSTH, timewindow, ExpType)

% written by MD
% function to bin PSTH
% Inputs:
% PSTH4D or 5D for 16 odors or conc exp
% timewindow = time in ms over which to bin
% ExpType = "Conc" if concentration exp or "Id" if 16 odors exp
% Smoothened PSTH with the following dimensions:
% Neurons x Odors X Concentration x timebin x Repeats if input was PSTH 5D
% Neurons x Stimuli x timebin x Repeats if input was PSTH 4D

Nneurons= size(PSTH,1);
clusterNum = 1:Nneurons;
%firstbin = -10;
%lastbin = 10;
t_wid = timewindow; %%specified by user in ms
%t_bin = t_wid./1000; % convert to sec
%edges = firstbin:t_bin:lastbin;
edges = 0:t_wid:20000;
x = 1:20000;
[~,~,idx] = histcounts(x, edges);

if ndims(PSTH) == 5 && ExpType == "Conc" % conc exp
    NTrials = 100;
    NOdors = 20;
    %timepoints = 20000;
elseif ndims(PSTH) == 4 && ExpType == "Conc" % conc exp and 4d PSTH matrix
    NTrials = 100;
    NOdors = 20;
elseif ndims(PSTH) == 4 && ExpType == "Id"  % 16 odors exp
    NTrials = 80;
    NOdors = 16;
    %timepoints = 20000;%in previous version of this exp it was shorter
end


%%
switch ndims(PSTH)
    %%
    case 5 % Compute smooth 5-D PSTH matrix as Neurons X Nber of Odors X Concentration X time (in ms) X Repeats
        
        NId = 5; % nber of odor id
        NConc = 4; % nber of od conc
        NRep = 5; % nber of repeats
        
        for c = 1:length(clusterNum) % loop through each unit
            clusterIdx = clusterNum(c);
            for x = 1:NId %for each odor identity
                for y = 1:NConc % for each concentration
                    for j = 1:5 % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                        tempPSTH = squeeze(PSTH(clusterIdx,x,y,:,j)); % for one cluster, all times for one type
                        %[myspikecount,idx] = histcounts(tempPSTH, edges);
                        %[~,~,idx] = histcounts(tempPSTH, edges);
                        %FR_mean_in_bin = accumarray(idx(2:101),myspikecount(:),[],@mean);
                        FR_mean_in_bin = accumarray(idx(:),tempPSTH(:),[],@mean);
                        binnedPSTH(clusterIdx,x,y,:,j) = FR_mean_in_bin;
                    end
                end
            end
        end
        %%
    case 4 % Compute 4-D PSTH as Neurons X Nber of Odors X time (in ms) X Repeats
        
        NRep = 5; % nber of repeats
        
        for c = 1:length(clusterNum) % loop through each unit
            clusterIdx = clusterNum(c);
            for x = 1:NOdors %for each odor identity
                for j = 1:5 % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                    tempPSTH = squeeze(PSTH(clusterIdx,x,:,j)); % for one cluster, all times for one type
                    FR_mean_in_bin = accumarray(idx(:),tempPSTH(:),[],@mean);
                    binnedPSTH(clusterIdx,x,:,j) = FR_mean_in_bin;
                end
            end
        end
        
end
end