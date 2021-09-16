function []= plotZScoredOdorConcMatrix(smoothPSTH)

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
tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]
firingRatesAverage = nanmean(x,5);

mu = squeeze(mean(firingRatesAverage(:,:,:,6000:10000),[2 3 4]));
sigma = squeeze(std(firingRatesAverage(:,:,:,6000:10000),0,[2 3 4]));
z_score = (firingRatesAverage(:,:,:,:) - mu) ./ sigma;

firstbin = -10;
lastbin = 10;

for whatsmell = 1:5
    figure(whatsmell)
    
    for whatconc = 1:4
        FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
        subplot(1,4,whatconc)
        X = squeeze(z_score(:,whatsmell,whatconc,:));
        h = imagesc(X);
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([-4  4]);
        xline(0,'k--');
        xline(4,'k--');
    end
    colormap(redblue);
    title( ['Odor',  num2str(whatsmell)])
    colorbar ('Position', [0.8 0.1 0.01 0.8]);
    
end

end