%PSTH_matrix_100clusters

%run prep_dPCA_matrix before 

goodcluster = allclusters;
NTrials = 100;
timepoints = 20000; 
Nneurons= length(goodcluster); % 97
step = 0.001; %timebin for PSTH should be 1 ms 

PSTH = zeros(Nneurons, NTrials,timepoints); %initialize 

%% Compute 3-D PSTH as Neurons X trials X time (in ms)

for whichcluster = 1:length(goodcluster)
    clear spikecount
    
    for mytrial = 1:NTrials
        firstbin = -10;
        lastbin = 10;
        tbins = firstbin:step:lastbin;
        SpikeTimes = goodcluster(whichcluster).spikes{1, mytrial};
        myspikecount = histcounts(SpikeTimes, tbins);
        PSTH(whichcluster,mytrial,:) = [myspikecount];
    end

end



%% Smooth PSTH 

taxis = -500:500;  % make a time axis of 1000 ms
t_wid = 100;  % width of kernel (in ms)
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

clusterNum = 1:Nneurons;
times = NTrials;

for c = 1:length(clusterNum) % loop through each unit
    clusterIdx = clusterNum(c);
    for t = 1:times % loop through trials
        tempPSTH = squeeze(PSTH(clusterIdx,t,:)); % for one cluster, all times for one trial 
        % squeeze enables having just one dimension

        %PSTH is where I store spike times (1 for spikes and 0 elsewhere)
        %PSTH is Neurons X trials X time (in ms)

        zs = conv(tempPSTH,gauss_kernel,'same');
        FR(clusterIdx,t,:) = zs*1000; %converting firing rate to Hz.

        % FR contains smoothed firing rates;

    end
end


%% Compute smooth 5-D PSTH matrix as Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 

NOdors = 20; %  for id/conc exp
NId = 5; % nber of odor id 
NConc = 4; % nber of od conc
NRep = 5; % nber of repeats 
PSTH5D = zeros(Nneurons, NId, NConc ,timepoints, NRep); %initialize 

Trial.Id(1,:) = [1 6 11 16]; %Odor 1
Trial.Id(2,:) = [2 7 12 17]; %Odor 2
Trial.Id(3,:) = [3 8 13 18]; %Odor 3
Trial.Id(4,:) = [4 9 14 19]; %Odor 4
Trial.Id(5,:) = [5 10 15 20]; % Odor 5
% 
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
                              % basically trial number of each rep for one
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
        
       % for j = 1:numel(reps) %repeat number % have to bypass that when >5
       % but will have to change for smaller 
       for j = 1:5
              
              % use these lines for unsmoothen data
%             firstbin = -10;
%             lastbin = 10;
%             tbins = firstbin:step:lastbin;
%             SpikeTimes = goodcluster(whichcluster).spikes{1, reps(j)};
%             myspikecount = histcounts(SpikeTimes, tbins);
%             PSTH5D(whichcluster,x,y,:,j) = [myspikecount];

              PSTH5D(whichcluster,x,y,:,j) = FR(whichcluster,reps(j),:);

        end

     end 
    
end


%% for checking purpose        

%using both matrices, cluster 1 trial 1 (in this case od 1, conc 1) all
%times
%figure(1); plot(squeeze(FR(1,1,:)))
%figure(2); plot(squeeze(PSTH5D(1,1,2,:,1)))

%using both matrices, cluster 23 trial 23 (in this case od 2, conc 3) all
%times
%figure(1); plot(squeeze(FR(23,23,:)))
%figure(2); plot(squeeze(PSTH5D(23,2,3,:,2)))

% cluster 1 all trials all times -> mean dim 2 to get all trials average
%figure(1); plot(squeeze(mean(FR(1,:,:),2))) 

% cluster 6 all trials all times -> mean dim 2 to get all trials average 
%figure(1); plot(squeeze(mean(FR(6,:,:),2))) % mean dim 2 

% cluster 6 odor 3 conc 1 all repeats -> mean dim 5 to get all repeats
% average of this stim type
%figure(2); plot(squeeze(mean(PSTH5D(6,3,1,:,:),5))) % mean dim 5

% cluster 6 all trials all times -> mean dim 2 to get all trials average
%figure(1); plot(squeeze(mean(FR(1,:,:),2)))