%% if PSTH5D is input

firstbin = -10;
lastbin = 10;

firingRates =  PSTH5D; %Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 
Nneurons= size(PSTH5D,1);
firingRatesAverage = nanmean(firingRates,5); % 4 DIM

for whatsmell =1:5
%whatsmell = 3;

FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:)); % For Odor 1 test


%myNeuron = [1 2];

%% neurons NOT organized by strength of response
figure(whatsmell) % NORMALIZED
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
        %normalizing 
        %airFiringRateAverage = squeeze(nanmean(FRAverageConc(myneuron,whatconc,1:4000),3));
        airFiringRateAverage = squeeze(nanmean(FRAverageConc(:,:,1:10000),[2 3]));
        Y = airFiringRateAverage;
        %Y = round(airFiringRateAverage);
        %X = ((FRAverageConc(myneuron,whatconc,:)) - Y) ./ Y
        Z = ((FRAverageConc(:,whatconc,:)) - Y) ./ Y;
        Z(isinf(Z)) = 0; %to avoid having inf values when getting a zero
        minvalueZ = min(Z,[],'all'); 
        maxvalueZ = max(Z,[],'all');
        allminZ(whatconc) = minvalueZ;
        allmaxZ(whatconc) = maxvalueZ;
        h = imagesc(squeeze(Z));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        %caxis([min(X,[],'all')  max(X,[],'all') ])
        caxis([-2  5]);
%        colormap(redblue);
        colormap(mymap);
%         xline(0,'w--');
%         xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear Z; clear Y;
end 

%end

% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])


%colormap(jet);

%% mean_plot

for whatconc = 1:4
        airmean = squeeze(nanmean(FRAverageConc(:,:,6000:10000),[1 2 3]));
        odormean = squeeze(mean(FRAverageConc(:,whatconc,10000:14000),[1 3]));
        stdeviation = squeeze(std(FRAverageConc(:,whatconc,10000:14000),0,[1 3]));
        normalized_mean = (odormean - airmean) / airmean;
        all_normalized_mean (whatconc,whatsmell) = normalized_mean;
end

X = [-4 -3 -2 -1];
Y = all_normalized_mean(:,whatsmell);
Y2 = mean(all_normalized_mean,2);

figure(6)
subplot(3,2,whatsmell)
plot (X,Y, '-b', 'LineWidth', 1.5)
xlim([-5 0])
ylim([min(all_normalized_mean(:,whatsmell))-0.05 max(all_normalized_mean(:,whatsmell)+0.05)])
xticks(X)
yline(0,'k--');

subplot(3,2,6)
plot (X,Y2, '-b', 'LineWidth', 1.5)
xlim([-5 0])
ylim([min(Y2)-0.05 max(Y2)+0.05])
xticks(X)
yline(0,'k--');

figure (12)
plot (X,Y, '-b', 'LineWidth', 1.5); hold on
xlim([-5 0])
xticks(X)
yline(0,'k--');

end

%% neurons ordered by strength of response

whatod = 1 ;%arbritarily chosing od 3 as ref for ordering
FRAverageConcOd3 = squeeze(firingRatesAverage(:,whatod,:,:)); % first attempt chosing od 3 as ref for ordering
airFiringRateAverage = squeeze(nanmean(FRAverageConcOd3(:,:,1:10000),[2 3]));
Normalized = ((FRAverageConcOd3(:,:,:)) - airFiringRateAverage) ./ airFiringRateAverage;
NormOdorfiringRatesAverage = squeeze(nanmean(Normalized(:,:,10000:14000),[2 3]));
[~, data_window_avg_idx]  = sort(NormOdorfiringRatesAverage,'descend'); 

for whatsmell =1:5
%whatsmell = 3;
FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
% % Get average response per neuron durinng stim
% 
% odorfiringRatesAverage = squeeze(nanmean(firingRatesAverage(:,:,:,10000:14000),[2 3 4]));
% 
%                               
% % Sorting of values from the window of interest
% 
% [~, data_window_avg_idx]  = sort(odorfiringRatesAverage,'descend'); 

% 
% % Sorting of ROI in original matrix (data) and average ROI response matrix (data_ROI_mean)
% 
% firingRatesAverage_ordered = firingRatesAverage(data_window_avg_idx,:,:,:); % data sorting
% FRAverageConc_ordered = squeeze(firingRatesAverage_ordered(:,whatsmell,:,:)); % For Odor 1 test


%myNeuron = [1 2];

figure (6+whatsmell) % NORMALIZED
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
        %normalizing 
        %airFiringRateAverage = squeeze(nanmean(FRAverageConc(myneuron,whatconc,1:4000),3));
        %airFiringRateAverage_ordered = squeeze(nanmean(FRAverageConc_ordered(:,:,1:10000),[2 3]));
        airFiringRateAverage = squeeze(nanmean(FRAverageConc(:,:,1:10000),[2 3]));
        %Y = airFiringRateAverage_ordered;
        Y = airFiringRateAverage;
        %Y = round(airFiringRateAverage);
        %X = ((FRAverageConc(myneuron,whatconc,:)) - Y) ./ Y
        %Z = ((FRAverageConc_ordered(:,whatconc,:)) - Y) ./ Y;
        Z = ((FRAverageConc(:,whatconc,:)) - Y) ./ Y;
        Z(isinf(Z)) = 0; %to avoid having inf values when getting a zero
%       
%         odorfiringRatesAverage = squeeze(nanmean(Z(:,:,10000:14000),[2 3]));
%         [~, data_window_avg_idx]  = sort(odorfiringRatesAverage,'descend'); 
        Z_ordered = Z(data_window_avg_idx,:,:); % data sorting
        
        %h = imagesc(squeeze(Z));
        h = imagesc(squeeze(Z_ordered));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        %caxis([min(X,[],'all')  max(X,[],'all') ])
        caxis([-2  5]);
%        colormap(redblue);
        colormap(mymap);
%         xline(0,'w--');
%         xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear Z; clear Y; clear Z_ordered;
end
end