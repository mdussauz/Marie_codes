function []= plotOdorConcActivityMatrix(smoothPSTH)

%written by MD
% Reponse Matrix for each concentration and each odor
%input PSTH 4D from conc experiments
% Smoothened PSTH with the following dimensions:
% Neurons x Stimuli x time x Repeats 
%plot odor-evoked activity averaged across repeats
%of all neurons for each conc and each odor
%-------------------------------------------------
tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc rep]
firingRatesAverage = nanmean(x,5);

firstbin = -10;
lastbin = 10;

for whatsmell = 1:5
    figure(whatsmell)
      
    for whatconc = 1:4
        subplot(1,5,whatconc)
        X = squeeze(firingRatesAverage(:,whatsmell,whatconc,:));
        h = imagesc(X);
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10]);
        caxis([0  70]);  %probably need to change range of legend and need to add legend at different spot
        xline(0,'k--');
        xline(4,'k--');
        clear X; clear Y;
    end
    colormap(parula);
    title( ['Odor',  num2str(whatsmell)])
    colorbar ('Position', [0.8 0.1 0.01 0.8]);
    
end

end