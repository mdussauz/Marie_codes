function[f] = correlationHeatMap(aircurves, odorcurves)
% Written by JW
% Plots the heatmaps between phase curves of air-odor and air-air
% Input: 
%   -Air phase curves = numStim x numBins
%   -Odor phase curves = numStim x numBins
% Output: figure
numStim = size(aircurves, 1);
[gridX, gridY] = getGrid(numStim);

f = figure();
colormap autumn
odorcorrelations = zeros(numStim, 1, 'double');
aircorrelations = zeros(numStim, 1, 'double');
for stim = 1:numStim
    odorcorr = corrcoef(aircurves(stim, :), odorcurves(stim, :));
    odorcorrelations(stim) = odorcorr(1, 2);
    aircorr = corrcoef(aircurves(stim, :), mean(odorcurves(:, :), 1));
    aircorrelations(stim) = aircorr(1, 2);
end
odorcorrelations = reshape(odorcorrelations, gridX, gridY);
aircorrelations = reshape(aircorrelations, gridX, gridY);
f = figure();
subplot(1, 2, 1)
imagesc(odorcorrelations, [-1 1]);
title("Air-Odor");
xlabel("Concentration");
ylabel("Identity");
subplot(1, 2, 2)
imagesc(aircorrelations, [-1 1]);
title("Air-Air");
xlabel("Concentration");
ylabel("Identity");
colorbar
sgtitle("Pearson Correlation between Phase Curves");
end