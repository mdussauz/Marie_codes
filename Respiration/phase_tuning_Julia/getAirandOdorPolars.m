function[airpolars, odorpolars] = getAirandOdorPolars(cluster, ind, phaseStarts)
% Written by JW
% Gets all the respiration phase polars for spikes for a particular
% stimulus for air and odor periods
% Input:
% -cluster: cluster from goodspikes
% -the stimulus index to plot
% -phaseStarts: an array of starts of respiration phases
spikes = cluster.spikes(ind, :);
odorspikes = cluster.odorspikes(ind, :);
spikes = vertcat(spikes{:});
odorspikes = vertcat(odorspikes{:});

%get spikes before
before = spikes(odorspikes < 0);
odor = spikes(odorspikes > 0.3 & odorspikes < 4);

airpolars = getAllPolars(before, phaseStarts);
odorpolars = getAllPolars(odor, phaseStarts);

end