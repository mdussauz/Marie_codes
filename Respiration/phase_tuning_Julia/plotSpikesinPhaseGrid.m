function[f] = plotSpikesinPhaseGrid(cluster, phaseStarts)
% plot spikes in phase for a cluster in a grid
% Each plot looks like : x-axis = time aligned to odor start, y-axis =
% respiration phase
% Input:
%   cluster
%   phaseStarts = a list of respiration phase start times
numStim = size(cluster.odorspikes, 1);
[gridX, gridY] = getGrid(numStim);
allspikes = cluster.spikes;
allodorspikes = cluster.odorspikes;
f = figure();
plot = 1;
allspikes = reshape(allspikes, gridX, gridY, 5);
allodorspikes = reshape(allodorspikes, gridX, gridY, 5);
for odorInd = 1:gridX
    for concInd = 1:gridY
        subplot(gridX, gridY, plot);
        for trialInd = 1:5
            spikes = allspikes{odorInd, concInd, trialInd};
            odorspikes = allodorspikes{odorInd, concInd, trialInd};
            y = zeros(size(spikes));
            for ind = 1:length(spikes)
                y(ind) = getPolar(spikes(ind), phaseStarts);
            end
            scatter(odorspikes, y, 0.2)
            xline(0)
            xline(4)
            set(gca, 'XLim', [-10, 10]);
            hold on
        end
        hold off
        plot = plot + 1;
    end
end
sgtitle("Raster plots");
end
