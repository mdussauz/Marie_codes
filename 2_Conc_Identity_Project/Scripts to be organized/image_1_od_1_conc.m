

%% if PSTH5D is input

firstbin = -10;
lastbin = 10;

firingRates =  PSTH5D; %Neurons X Nber of Odors X Concentration X time (in ms) X Repeats 
Nneurons= size(PSTH5D,1);

% NORMALIZING
% airFiringRateAverage = nanmean(firingRates(:,:,:,1:10000,:),4);
% normalizedFiringRates = (firingRates - airFiringRateAverage)./airFiringRateAverage;
% normalizedFiringRatesAverage = nanmean(normalizedFiringRates,5);
% normalizedFRAverageConc = squeeze(mean(normalizedFiringRatesAverage, 2));
% normalizedFRAverageOd = squeeze(mean(normalizedFiringRatesAverage, 3));



%WITHOUT NORMALIZING
whatsmell = 5;
firingRatesAverage = nanmean(firingRates,5);
FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:)); % For Odor 1 test
% FRAverageConc = squeeze(mean(firingRatesAverage, 2));
% FRAverageOd = squeeze(mean(firingRatesAverage, 3));


%% New Plot 
%myNeuron = [1 2];
figure(1) 
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
     %   X = FRAverageConc(myneuron,whatconc,:);
        X = FRAverageConc(:,whatconc,:);
        minvalue = min(X,[],'all'); 
        maxvalue = max(X,[],'all');
        allmin(whatconc) = minvalue;
        allmax(whatconc) = maxvalue;
        h = imagesc(squeeze(X));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10]);
        caxis([0  70]);
        %caxis([min(X,[],'all')  max(X,[],'all') ])
%         xline(0,'w--');
%         xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear X; clear Y;
end 

%end

% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])


colormap(jet);
%colormap(gray);
%caxis([min(allmin)  max( allmax) ]);
%subplot(1,5,5)
%colorbar

figure(2) % NORMALIZED
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
        caxis([-5  5]);
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
colormap(redblue);
%caxis([min(allminZ)  max( allmaxZ) ]);

%subplot(1,5,5)
%colorbar

%%
figure(3) % Z SCORE
%for myNeuron = 1:Nneurons %for each neuron
for whatconc = 1:4
        subplot(1,4,whatconc)
        mu = squeeze(mean(FRAverageConc(:,:,1:10000),[2 3]));
        sigma = squeeze(std(FRAverageConc(:,:,1:10000),0,[2 3]));
        z_score = (FRAverageConc(:,whatconc,:) - mu) ./ sigma;
        h = imagesc(squeeze(z_score));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        %caxis([min(X,[],'all')  max(X,[],'all') ])
        caxis([-10  10]);
%         xline(0,'w--');
%         xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear mu; clear sigma; clear z_score;
end 

%end

% title( ['Response matrix for all trials of cluster #', ...
%     num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step)])


colormap(redblue);