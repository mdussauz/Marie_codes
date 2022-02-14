%% Smooth PSTH 

taxis = -500:500;  % make a time axis of 1000 ms
t_wid = 200;  % width of kernel (in ms)
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
        
        for j = 1:numel(reps) %repeat number 
              
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
%%

firstbin = -10;
lastbin = 10;

Nneurons= length(goodcluster);

firingRates =  PSTH5D; %Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 

% NORMALIZING
% airFiringRateAverage = nanmean(firingRates(:,:,:,1:10000,:),4);
% normalizedFiringRates = (firingRates - airFiringRateAverage)./airFiringRateAverage;
% normalizedFiringRatesAverage = nanmean(normalizedFiringRates,5);
% normalizedFRAverageConc = squeeze(mean(normalizedFiringRatesAverage, 2));
% normalizedFRAverageOd = squeeze(mean(normalizedFiringRatesAverage, 3));



%WITHOUT NORMALIZING
firingRatesAverage = nanmean(firingRates,5);
FRAverageConc = squeeze(mean(firingRatesAverage, 2));
FRAverageOd = squeeze(mean(firingRatesAverage, 3));






%% Plot 

figure(1) 
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
        X = FRAverageConc(:,whatconc,:);
        h = imagesc(squeeze(X))
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([min(X,[],'all')  max(X,[],'all') ])
        xline(0,'w--');
        xline(4,'w--');
        clear X; clear Y;
end 

%end

% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])


colormap(jet)
%subplot(1,5,5)
%colorbar

figure(2) % NORMALIZED
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
        %normalizing 
        airFiringRateAverage = nanmean(FRAverageConc(:,whatconc,1:4000),3);
        %Y = airFiringRateAverage;
        Y = round(airFiringRateAverage);
        X = ((FRAverageConc(:,whatconc,:)) - Y) ./ Y;
        X(isinf(X)) = 0; 
        h = imagesc(squeeze(X))
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([min(X,[],'all')  max(X,[],'all') ])
        xline(0,'w--');
        xline(4,'w--');
        clear X; clear Y;
end 

%end

% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])


colormap(jet)
%subplot(1,5,5)
%colorbar

figure(3) 
%for myNeuron = 1:Nneurons %for each neuron
for whatod = 1:5
        subplot(1,5,whatod)
        X = FRAverageOd(:,whatod,:);
        g = imagesc(squeeze(X))
        set(g, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([min(X,[],'all')  max(X,[],'all') ])
        xline(0,'w--');
        xline(4,'w--');
        clear X; clear Y;
end 

%end
        
% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])

colormap(jet)
%colorbar

figure(4)% NORMALIZED
%for myNeuron = 1:Nneurons %for each neuron
for whatod = 1:5
        subplot(1,5,whatod)
        %normalizing 
        airFiringRateAverage = nanmean(FRAverageOd(:,whatod,1:10000),3);
        Y = round(airFiringRateAverage);
        X = ((FRAverageOd(:,whatod,:)) - Y) ./ Y;
        X(isinf(X)) = 0; 
        g = imagesc(squeeze(X))
        set(g, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([min(X,[],'all')  max(X,[],'all') ])
        xline(0,'w--');
        xline(4,'w--');
        clear X; clear Y;
end 

%end
        
% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])

colormap(jet)
%colorbar


