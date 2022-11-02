function []= plotZScoredOdorConcMatrix_new(smoothPSTH)

%written by MD
% Z SCORE response matrix for each concentration and each odor
%input PSTH 4D from conc experiments
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
firingRatesAverage = nanmean(smoothPSTH,5);

mu = squeeze(mean(firingRatesAverage(:,:,:,17000:21000),[2 4])); %to be checked
sigma = squeeze(std(firingRatesAverage(:,:,:,17000:21000),0,[2 4])); %to be checked 
z_score = (firingRatesAverage - mu) ./ sigma;

firstbin = -21;
lastbin = 19;

for whatsmell = 1:5
    figure(whatsmell)
    
    for whatconc = 1:4
        subplot(1,4,whatconc)
        X = squeeze(z_score(:,whatsmell,whatconc,:));
        h = imagesc(X);
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([-4  4]);
        xline(0,'k--');
        xline(2,'k--');
    end
    colormap(redblue);
    title( ['Odor',  num2str(whatsmell)])
    colorbar ('Position', [0.95 0.1 0.01 0.8]);
    
end

end