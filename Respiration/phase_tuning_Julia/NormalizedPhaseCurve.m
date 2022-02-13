function[hist, edges] = NormalizedPhaseCurve(polars, nbins)
% Gets a histogram of a normalized phase curve
% Input: 
% - polars = a list of spikes as polars
% - nbins = number of bins
    [hist, edges] = histcounts(polars, nbins, 'Normalization', 'probability');
end
