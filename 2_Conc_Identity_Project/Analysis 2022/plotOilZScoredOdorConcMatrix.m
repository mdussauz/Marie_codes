function [mu,sigma, firingRatesAverage]= plotOilZScoredOdorConcMatrix(smoothPSTH)

%written by MD
% Z SCORE response matrix for each concentration and each odor
% input PSTH 5D from conc experiments
% Smoothened PSTH with the following dimensions:
%  Neurons x identity x concentration x time x Repeats 
% plot z scored odor-evoked activity averaged across repeats
% of all neurons for each conc and each odor

% According to Schnoonover et al 2020, z scoring the smoothed PSTH by 
% substracting for each single unit its spontaneous baseline firing rate 
% across all trials (4s before stimulus onset) and divided by the 
% standard deviation of firing rates during the baseline epoch.
% Hence response maps indicates changes in FR in units of std of 
% spontaneous activity. 
%-------------------------------------------------
global prestim
global odorstim
global poststim

firstbin = -prestim/1000;
lastbin = (odorstim+poststim)/1000;

firingRatesAverage = mean(smoothPSTH,5); %dim: neuron x odor x concentration x time
FRAverageOil = squeeze(smoothPSTH(:,5,:,:,:));

mu = squeeze(mean(FRAverageOil(:,:,prestim:(prestim+odorstim),:),[3 4])); 
sigma = squeeze(std(FRAverageOil(:,:,prestim:(prestim+odorstim),:),0,[3 4]));  


for whatsmell = 1:5
    figure(whatsmell)
    
    for whatconc = 1:4
        subplot(1,4,whatconc)
        xx = squeeze(firingRatesAverage(:,whatsmell,whatconc,prestim:(prestim+odorstim)));
        z_score = ( xx - mu(:,whatconc)) ./ sigma(:,whatconc);
        z_score(isinf(z_score==inf))=0;
        z_score(isnan(z_score))=0;
        X = z_score;
        h = imagesc(X);
        set(h, 'XData', [0, odorstim/1000]);
        xlim([0, odorstim/1000])
        caxis([-2 2]);
    end
    colormap(redblue);
    title( ['Odor',  num2str(whatsmell)])
    colorbar ('Position', [0.95 0.1 0.01 0.8]);
    
end

end