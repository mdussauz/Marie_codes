function [f] = plotPSTHsinGrid(cluster,t_wid)
%%
%written by MD
%plot PSTHs averaged across repeats in a grid in a grid 
% for a subplot 
% - y axis = Firing Rate
% - x axis = time in sec
% for the grid in the case of concentration exp
% - y axis = different odors
% - x axis = concentration 
% subplot(1) = odor 1, conc 10^-4
% subplot(20) = odor 5, conc 10^-1

%%
%allodorspikes = cluster.odorspikes;
f = figure();
plots = 1;
numStim = size(cluster.odorspikes, 1);

%%
[smoothPSTH,tbins] = ComputePSTH(cluster,t_wid);
PSTH = smoothPSTH;
%%
function[gridX, gridY] = getGrid(numStim)
% Gives gridX and gridY if 20 stimuli (concentration) or 16 stimuli
% (identity)
    if numStim == 20
        gridX = 5;
        gridY = 4;
    else
        gridX = 4;
        gridY = 4;
    end
end

%%
[gridX, gridY] = getGrid(numStim);
%% Plot
%PSTH as Nber of Odors X Concentration X time (in ms) X Repeats 
    for odorInd = 1:gridX
        for concInd = 1:gridY
            subplot(gridX, gridY, plots);
            RepAveragePSTH = squeeze(mean(PSTH(odorInd,concInd,:,:),4));
                
            plot(tbins(1:end-1), RepAveragePSTH)
            xline(0);
            xline(4);
            set(gca, 'XLim', [-10, 10]);
            hold off
            plots = plots + 1; 
        end
 
    end
    
sgtitle("PSTH plots");
end