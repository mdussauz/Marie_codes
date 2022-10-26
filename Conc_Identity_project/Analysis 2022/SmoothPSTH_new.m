function [smoothPSTH] = SmoothPSTH_new(PSTH, timewindow, ExpType)

% written by MD
% function to smooth PSTH
% Inputs:
% PSTH4D or 5D for 16 odors or conc exp
% timewindow = time in ms over which to smooth 
% ExpType = "Conc" if concentration exp or "Id" if 16 odors exp
% Output:
% Smoothened PSTH with the following dimensions:
% Neurons x Odors X Concentration x time x Repeats if input was PSTH 5D
% Neurons x Stimuli x time x Repeats if input was PSTH 4D

Nneurons= size(PSTH,1);
clusterNum = 1:Nneurons;
t_wid = timewindow;  % width of kernel (in ms)
taxis = -(t_wid*5):(t_wid*5);  % make a time axis of 1000 ms for a window of 100 ms
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

if ndims(PSTH) == 5 && ExpType == "Conc" % conc exp
    NTrials = 140;
    NOdors = 20;
	%timepoints = 20000; 
elseif ndims(PSTH) == 4 && ExpType == "Conc" % conc exp and 4d PSTH matrix
    NTrials = 140;
    NOdors = 20;
elseif ndims(PSTH) == 4 && ExpType == "Id"  % 16 odors exp
    NTrials = 112;
    NOdors = 16;
    %timepoints = 20000;%in previous version of this exp it was shorter    
end
   

%%
switch ndims(PSTH)
%%
    case 5 % Compute smooth 5-D PSTH matrix as Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 

    NId = size(PSTH,2); % nber of odor id 
    NConc = size(PSTH,3); % nber of od conc
    NRep = size(PSTH,5); % nber of repeats 
    
    for c = 1:length(clusterNum) % loop through each unit
        clusterIdx = clusterNum(c);
        for x = 1:NId %for each odor identity
            for y = 1:NConc % for each concentration
                for j = 1:NRep % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                    tempPSTH = squeeze(PSTH(clusterIdx,x,y,:,j)); % for one cluster, all times for one type 
                    zs = conv(tempPSTH,gauss_kernel,'same');
                    smoothPSTH(clusterIdx,x,y,:,j) = zs*1000; %converting firing rate to Hz.
                end
            end
        end 
    end
%%
    case 4 % Compute 4-D PSTH as Neurons X Nber of Stimuli X time (in ms) X Repeats 

    NRep = size(PSTH,5); % nber of repeats  % nber of repeats 
    
    for c = 1:length(clusterNum) % loop through each unit
        clusterIdx = clusterNum(c);
        for x = 1:NOdors %for each odor stimulus
            for j = 1:NRep % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                tempPSTH = squeeze(PSTH(clusterIdx,x,:,j)); % for one cluster, all times for one type 
                zs = conv(tempPSTH,gauss_kernel,'same');
                smoothPSTH(clusterIdx,x,:,j) = zs*1000; %converting firing rate to Hz.
            end
        end 
    end

end 
end