function[beforecurve, odorcurve]= PhaseCurves(cluster, ind, phaseStarts,nbins)
% Plots phase air and odor phase curves in one plot
% Input:
%   - cluster = the cluster
%   - ind = the stimulus index
%   - phaseStarts = a list of respiration phase starts
%   - nbins = number of bins

[beforepolars, odorpolars] = getAirandOdorPolars(cluster, ind, phaseStarts);

[beforecurve, beforeedges] = NormalizedPhaseCurve(beforepolars, nbins);
plot(beforecurve);
hold on
[odorcurve, odoredges] = NormalizedPhaseCurve(odorpolars, nbins);
plot(odorcurve);
hold off
end

