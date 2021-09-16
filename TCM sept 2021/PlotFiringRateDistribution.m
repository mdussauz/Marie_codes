function [] = PlotFiringRateDistribution(smoothPSTH)

%written by MD
%visualize distribution of mean FR during baseline epoch 

tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]
firingRatesAverage = nanmean(x,5);
baselineFR = squeeze(mean(firingRatesAverage(:,:,:,6000:10000),[2 3 4]));

binsize = size(unique(baselineFR));
figure(1); histogram(baselineFR, binsize)
figure (2); boxplot(baselineFR)

end 