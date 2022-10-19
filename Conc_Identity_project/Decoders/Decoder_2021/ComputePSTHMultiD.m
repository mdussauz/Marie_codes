function [PSTH] = ComputePSTHMultiD(allclusters, ExpType,Dim)
%written by MD
% inputs:
% allclusters = structure array produced by dPCA prep code
% ExpType = "Conc" if concentration exp or "Id" if 16 odors exp
% Dim = 4 for PSTH4d or 5 for PSTH5d
% Outputs:
% PSTH with the following dimensions:
% Neurons x Odors X Concentration x time x Repeats - for concentration series exp only and if we want to separate
% identity and concentration in different dimensions 
% Neurons x Stimuli x time x Repeats -  for both conc series and 16 odors
% exp 

goodcluster = allclusters;
firstbin = -10;
lastbin = 10;
step = 0.001; %timebin for PSTH should be 1 ms 
tbins = firstbin:step:lastbin;
Nneurons= length(goodcluster); 

if ExpType == "Conc"
    %NTrials = 100;
    NOdors = 20;
	timepoints = 20000; 
elseif ExpType == "Id" 
    %NTrials = 80;
    NOdors = 16;
    timepoints = 20000;%in previous version of this exp it was shorter    
end
    
%%
switch Dim
%%    
    case 5 % Compute 5-D PSTH matrix as Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 

    NId = 5; % nber of odor id 
    NConc = 4; % nber of od conc
    NRep = 5; % nber of repeats 
    PSTH5D = zeros(Nneurons, NId, NConc ,timepoints, NRep); %initialize 

    Trial.Id(1,:) = [1 6 11 16]; %Odor 1
    Trial.Id(2,:) = [2 7 12 17]; %Odor 2
    Trial.Id(3,:) = [3 8 13 18]; %Odor 3
    Trial.Id(4,:) = [4 9 14 19]; %Odor 4
    Trial.Id(5,:) = [5 10 15 20]; % Odor 5

    Trial.Conc(1,:) = [1 2 3 4 5];% concentration 10^-4 
    Trial.Conc(2,:) = [6 7 8 9 10];% concentration 10^-3 
    Trial.Conc(3,:) = [11 12 13 14 15];% concentration 10^-2 
    Trial.Conc(4,:) = [16 17 18 19 20];% concentration 10^-1 

    for whichcluster = 1:length(goodcluster)
    stim = goodcluster(whichcluster).stimulus;
    clear spikecount
    
        for i = 1:NOdors %for each odor % 20 for id/conc exp
            clear spikecount 
            reps = find(stim==i); % get indices of all repeats for each odor  
                              % ~trial number of each rep for one
                              % odor
            if ismember(i , Trial.Id(1,:))
                x = 1;
            elseif ismember(i , Trial.Id(2,:))
                x = 2;
            elseif ismember(i , Trial.Id(3,:))
                x = 3;
            elseif ismember(i , Trial.Id(4,:))    
                x = 4;
            elseif ismember(i , Trial.Id(5,:))    
                x = 5;
            end 
        
            if ismember(i ,Trial.Conc(1,:)) 
                y = 1;
            elseif ismember(i ,Trial.Conc(2,:)) 
                y = 2;
            elseif ismember(i ,Trial.Conc(3,:)) 
                y = 3;
            elseif ismember(i ,Trial.Conc(4,:)) 
                y = 4;
            end 

           for j = 1:5 % for j = 1:numel(reps) %repeat number % have to bypass that when >5
        SpikeTimes = goodcluster(whichcluster).spikes{1, reps(j)};
        myspikecount = histcounts(SpikeTimes, tbins); %partitions the spiketimes into bins, and returns the count in each bin
        PSTH5D(whichcluster,x,y,:,j) = [myspikecount];
           end
        end 
    end
    PSTH = PSTH5D;
%%
    case 4 % Compute 4-D PSTH as Neurons X Nber of Odors X time (in ms) X Repeats 

    NRep = 5; % nber of repeats 
    PSTH4D = zeros(Nneurons, NOdors ,timepoints, NRep); %initialize 


    for whichcluster = 1:length(goodcluster)
        stim = goodcluster(whichcluster).stimulus;
        clear spikecount
        for i = 1:NOdors %for each odor % 20 for id/conc exp or 16 for 16 od exp
            clear spikecount 
            reps = find(stim==i); % get indices of all repeats for each odor  
                              % ~trial number of each rep for one
                              % odor
            for j = 1:5
                SpikeTimes = goodcluster(whichcluster).spikes{1, reps(j)};
                myspikecount = histcounts(SpikeTimes, tbins);
                PSTH4D(whichcluster,i,:,j) = [myspikecount];
            end
        end 
    end
    PSTH = PSTH4D;
end 
    
end



