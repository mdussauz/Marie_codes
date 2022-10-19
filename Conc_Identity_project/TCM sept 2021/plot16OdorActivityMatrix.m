function []= plot16OdorActivityMatrix(smoothPSTH)

%written by MD
% Reponse Matrix for each concentration and each odor
%input PSTH 4D from conc experiments
% Smoothened PSTH with the following dimensions:
% Neurons x Stimuli x time x Repeats
%plot odor-evoked activity averaged across repeats
%of all neurons for each conc and each odor
%-------------------------------------------------
firingRatesAverage = nanmean(smoothPSTH,4);
firstbin = -10;
lastbin = 10;

figCount = 1;
subplotcount = 1;

for whatsmell = 1:16
    figure(figCount)
    subplot(1,4,subplotcount)
    
    X = squeeze(firingRatesAverage(:,whatsmell,:));
    h = imagesc(X);
    set(h, 'XData', [firstbin, lastbin]);
    xlim([-10, 10]);
    caxis([0  30]);  %probably need to change range of legend and need to add legend at different spot
    xline(0,'k--');
    xline(4,'k--');
    clear X; clear Y;
    
    colormap(parula);
    title( ['Odor',  num2str(whatsmell)])
    
    subplotcount = subplotcount +1;
    if subplotcount == 5
        subplotcount =1;
        figCount = figCount +1;
    end
    
    
end
colorbar ('Position', [0.8 0.1 0.01 0.8]);
end