function[f] = FiringRateGrid(cluster)
% Written by JW
% Plots the firing rate as bar graph of air and odor in a grid
% Input: 
%   - cluster: the cluster from goodcluster structure
% Output: figure
numStim = size(cluster.odorspikes, 1);
[gridX, gridY] = getGrid(numStim);
f = figure();
indConverter = reshape(1:numStim, gridX, gridY);
indConverter = reshape(indConverter', numStim, 1);
for stim = 1:numStim
    
    subplot(gridX, gridY, stim)
    odorspikes = cluster.odorspikes(indConverter(stim), :);
    odorspikes = vertcat(odorspikes{:});
    air = sum(odorspikes < 0)/50;
    odor = sum(odorspikes > 0.3 & odorspikes < 4)/20;
    bar(0, air);
    hold on 
    bar(1,odor );
    
end
sgtitle("Average Firing Rate");

end