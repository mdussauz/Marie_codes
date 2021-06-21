function [smoothPSTH] = SmoothPSTH(PSTH, timewindow, ExpType)

% written by MD
% function to smooth PSTH
% PSTH4D or 5D for 16 odors or conc exp
% timewindow = time in ms over which to smooth 
% ExpType = "Conc" if concentration exp or "Id" if 16 odors exp

Nneurons= size(PSTH,1);
t_wid = timewindow;  % width of kernel (in ms)
taxis = -(t_wid*5):(t_wid*5);  % make a time axis of 1000 ms for a window of 100 ms
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);
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
                for j = 1:5 % for j = 1:numel(reps) %repeat number % have to bypass that when >5
                    tempPSTH = squeeze(PSTH(clusterIdx,x,y,:,j)); % for one cluster, all times for one type 
                    zs = conv(tempPSTH,gauss_kernel,'same');
                    smoothPSTH(clusterIdx,x,y,:,j) = zs*1000; %converting firing rate to Hz.
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
                zs = conv(tempPSTH,gauss_kernel,'same');
                smoothPSTH(clusterIdx,x,:,j) = zs*1000; %converting firing rate to Hz.
            end
        end 
    end

end 
end