function[polars] = getAllPolars(sp, phaseStarts)
% Transforms a list of spikes into a list of respiration phase polars
% Input:
% - sp = list of spikes
% - phaseStarts = a list of respiration phase starts
polars = zeros(length(sp), 1, 'double');
for ind = 1:length(sp)
    polars(ind) = getPolar(sp(ind), phaseStarts);
    
end

end