
function [normalizedPSTH] = NormalizedPSTH(PSTH, ExpType)

% written by MD
% function to nromalize PSTH by substracting baseline
% Inputs: 
% PSTH4D or 5D for 16 odors or conc exp
% ExpType = "Conc" if concentration exp or "Id" if 16 odors exp
% Output:
% Normalized PSTH with the following dimensions:
% Timebin x Neurons x Odors X Concentration x Repeats if input was PSTH 5D
% Timebin x Neurons x Stimuli X Repeats  if input was PSTH 5D


Nneurons= size(PSTH,1);
clusterNum = 1:Nneurons;

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
                    for j = 1:NRep % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                            tempPSTH = squeeze(PSTH(clusterIdx,x,y,:,j)); % for one cluster, all times for one type
                            airmean = squeeze(nanmean(tempPSTH(9000:10000)));
                            temp_normalized_PSTH = (tempPSTH - airmean); % / airmean;
                            clear airmean;
                            normalizedPSTH(:,clusterIdx,x,y,j) = temp_normalized_PSTH;
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
                for j = 1:NRep % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                        tempPSTH = squeeze(PSTH(clusterIdx,x,:,j)); % for one cluster, all times for one type
                        airmean = squeeze(nanmean(tempPSTH(9000:10000)));
                        temp_normalized_PSTH = (tempPSTH - airmean); % / airmean;
                        clear airmean;
                        normalizedPSTH(:,clusterIdx,x,j) = temp_normalized_PSTH;
                end
            end
            
        end
        
end
end