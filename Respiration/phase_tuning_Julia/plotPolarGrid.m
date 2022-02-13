function[f, deltaThetas] = plotPolarGrid(cluster, phaseStarts)
% plots the average air and odor vectors in a grid for a cluster
    numStim = size(cluster.odorspikes, 1);
    [gridX, gridY] = getGrid(numStim);
    f = figure();
    deltaThetas = zeros(numStim, 1, 'double');
    
    indConverter = reshape(1:numStim, gridX, gridY);
    indConverter = reshape(indConverter', numStim, 1);
    for stim = 1:numStim
        subplot(gridX, gridY, stim);
        [deltaTheta] = plotAirandOdorPolar(cluster, indConverter(stim), phaseStarts); 
        deltaThetas(indConverter(stim)) = deltaTheta;
    end
    sgtitle("Average vectors");
end

