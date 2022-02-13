function[f, allair, allodor] = phaseCurvesGrid(cluster, phaseStarts, nbins)
% Plot the air and odor phase curves for a cluster for all stimuli in a grid
% - Input:
%   - cluster
%   - phaseStarts = a list of respiration phase starts
%   - nbins = number of bins for the phase curves

numStim = size(cluster.odorspikes, 1);
[gridX, gridY] = getGrid(numStim);

allair = zeros(numStim, nbins, 'double');
allodor = zeros(numStim, nbins, 'double');
f = figure();
plotInd = 1;
indConverter = reshape(1:numStim, gridX, gridY);
indConverter = reshape(indConverter', numStim, 1);

for odorInd = 1:gridX
    for concInd = 1:gridY
        subplot(gridX, gridY, plotInd);
        [beforecurve, odorcurve] = PhaseCurves(cluster, indConverter(plotInd), phaseStarts, nbins);
        allair(indConverter(plotInd), :) = beforecurve;
        allodor(indConverter(plotInd), :) = odorcurve;
        plotInd = plotInd + 1;
        if odorInd == 1 && concInd == 1
            legend('before', 'odor');
        end
    end
    
end
sgtitle("Phase Curves");
end
