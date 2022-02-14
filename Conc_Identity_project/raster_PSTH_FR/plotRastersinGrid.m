function [f] = plotRastersinGrid(cluster,fast)
%%
%written by MD
%plot rasters in a grid 
% for a subplot 
% - y axis = 5 repeats, from bottom to top
% - x axis = time in sec
% for the grid in the case of concentration exp
% - y axis = different odors
% - x axis = concentration 
% subplot(1) = odor 1, conc 10^-4
% subplot(20) = odor 5, conc 10^-1

%%
allodorspikes = cluster.odorspikes;
f = figure();
plots = 1;
numStim = size(cluster.odorspikes, 1);
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
allodorspikes = reshape(allodorspikes, gridX, gridY, 5);

%%
switch fast 
    case 1 % option to plot spike as dots
    for odorInd = 1:gridX
        for concInd = 1:gridY
            subplot(gridX, gridY, plots);
            for trialInd = 1:5
                odorspikes = allodorspikes{odorInd, concInd, trialInd};
                y = zeros(size(odorspikes));
                y = y+trialInd/2;
                scatter(odorspikes,y, 0.2)
                xline(0);
                xline(4);
                set(gca, 'XLim', [-10, 10]);
                hold on
            end
            hold off
            plots = plots + 1;
        end
    end


case 0 % option to plot spike as lines
    for odorInd = 1:gridX
        for concInd = 1:gridY
            subplot(gridX, gridY, plots);
            for trialInd = 1:5
                odorspikes = allodorspikes{odorInd, concInd, trialInd};
                x = [1 1].*odorspikes; 
                y = ones(size(odorspikes));
                y = [0 1].*y; 
                y = y +trialInd-1;
     
                for spikes = 1: numel(odorspikes)
                    line(x(spikes,:),y(spikes,:))
                end
                xline(0); % starts of odor period
                xline(4); % end of odor period
                set(gca, 'XLim', [-10, 10]);
                set(gca, 'YLim', [0, 5]);
                hold on
                 
            end
            hold off
            plots = plots + 1;
        end
    end

end 
    
sgtitle("Raster plots");
end
    