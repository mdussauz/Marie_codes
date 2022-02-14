function [smoothPSTH, tbins] = ComputePSTH(cluster,t_wid)
%%
%written by MD
%function to compute PSTH and smoothen it with a gaussian filter
%%
allodorspikes = cluster.odorspikes;
numStim = size(cluster.odorspikes, 1);
repeats = 5;
NTrials = numStim*repeats;
timepoints = 20000;

firstbin = -10;
lastbin = 10;
step = 0.001; %timebin for PSTH should be 1 ms 
tbins = firstbin:step:lastbin;

%%
function[gridX, gridY] = getGrid(numStim)
% Gives gridX and gridY if 20 stimuli (concentration) or 16 stimuli
% (identity)
    if numStim == 20
        gridX = 5;
        gridY = 4;
    else
        gridX = 4;
        gridY = 4;
    end
end

%%
[gridX, gridY] = getGrid(numStim);
allodorspikes = reshape(allodorspikes, gridX, gridY, 5);

%%
PSTH = zeros(gridX, gridX,timepoints,repeats); %initialize 

%% Compute PSTH as Nber of Odors X Concentration X time (in ms) X Repeats 

    for odorInd = 1:gridX
        for concInd = 1:gridY
            for trialInd = 1:5
                odorspikes = allodorspikes{odorInd, concInd, trialInd};
                myspikecount = histcounts(odorspikes, tbins);
                PSTH(odorInd, concInd,:, trialInd) = [myspikecount];
            end

        end
    end


%% Smooth PSTH 

taxis = -500:500;  % make a time axis of 1000 ms
%t_wid = 100;  % width of kernel (in ms)
%t_wid = 200; % trying a different one
%t_wid = 25; % trying a different one
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);


for odorInd = 1:gridX
    for concInd = 1:gridY
        for trialInd = 1:5
            tempPSTH = squeeze(PSTH(odorInd, concInd,:, trialInd)); 
          % for one cluster, all times for one trial 
          % squeeze enables having just one dimension
            zs = conv(tempPSTH,gauss_kernel,'same');
            smoothPSTH(odorInd, concInd,:, trialInd) = zs*1000; 
            %converting firing rate to Hz.
        end
    end
end

   
end