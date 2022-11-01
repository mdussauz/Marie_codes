function []= plotZScored16OdorMatrix_new(smoothPSTH)

%written by MD
% Z SCORE response matrix for each concentration and each odor
%input PSTH 4D from 16 odors experiments
% Smoothened PSTH with the following dimensions:
% Neurons x Stimuli x time x Repeats
%plot z scored odor-evoked activity averaged across repeats
%of all neurons for each conc and each odor

% According to Schnoonover et al 2020, z scoring the smoothed PSTH by
% substracting for each single unit its spontaneous baseline firing rate
% across all trials (4s before stimulus onset) and divided by the
% standard deviation of firing rates during the baseline epoch.
% Hence response maps indicates changes in FR in units of std of
% spontaneous activity.
%-------------------------------------------------
firingRatesAverage = nanmean(smoothPSTH,4);

mu = squeeze(mean(firingRatesAverage(:,:,17000:21000),[2 3]));
sigma = squeeze(std(firingRatesAverage(:,:,17000:21000),0,[2 3]));
z_score = (firingRatesAverage(:,:,:) - mu) ./ sigma;

firstbin = -21;
lastbin = 19;

figCount = 1;
subplotcount = 1;

for whatsmell = 1:16
    figure(figCount)
    subplot(1,4,subplotcount)
    
    X = squeeze(z_score(:,whatsmell,:));
    h = imagesc(X);
    set(h, 'XData', [firstbin, lastbin]);
    xlim([firstbin, lastbin])
    caxis([-4  4]);
    xline(0,'k--');
    xline(2,'k--');
    
    colormap(redblue);
    title( ['Odor',  num2str(whatsmell)])
    
    subplotcount = subplotcount +1;
    if subplotcount == 5
        subplotcount =1;
        figCount = figCount +1;
    end
    
end

colorbar ('Position', [0.95 0.1 0.01 0.8]);
end