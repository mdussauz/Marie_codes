function[deltaTheta] = plotAirandOdorPolar(cluster, ind, phaseStarts)
% plots the average air and odor vector for a stimulus for a cluster
% Input:
%   - cluster
%   - ind : stimulus index
%   - phaseStarts: a list of respiration phase starts

[airpolars, odorpolars] = getAirandOdorPolars(cluster, ind, phaseStarts);
[airTheta, airRho] = averagePolar(airpolars);
polarplot([0 airTheta], [0 airRho], 'Color', [15 89 238]/256, 'LineWidth', 2);
hold on
[odorTheta, odorRho] = averagePolar(odorpolars);
polarplot([0 odorTheta], [0 odorRho], 'Color', [240 71 71]/256, 'LineWidth', 2);
hold off
deltaTheta = abs(airTheta - odorTheta);
Ax = gca; % current axes
Ax.ThetaGrid = 'off';
Ax.RGrid = 'off';
Ax.RTickLabel = []; 
Ax.ThetaTickLabel = [];
rlim([0 1]);
end
