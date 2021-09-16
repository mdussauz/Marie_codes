function []= plotRepExampleOdorConcMatrix(smoothPSTH, repeat)

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
ExampleRepActivity = squeeze(x(:,:,:,:,repeat));

firstbin = -10;
lastbin = 10;

for whatsmell = 1:5
    figure(whatsmell)
      
    for whatconc = 1:4
        subplot(1,4,whatconc)
        X = squeeze(ExampleRepActivity(:,whatsmell,whatconc,:));
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
    
end

end