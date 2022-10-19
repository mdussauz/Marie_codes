function [] = plotmeanresponse16od(smoothPSTH)

firingRatesAverage = nanmean(smoothPSTH,4);

mean_FR = squeeze(mean(firingRatesAverage(:,:,10000:14000), 3));
plot1 = subplot(1,2,1)
h1 = imagesc(mean_FR);
caxis([0  30]); %use 20 for APC
colormap(plot1,parula);
colorbar

mu = squeeze(mean(firingRatesAverage(:,:,6000:10000),[2 3]));
sigma = squeeze(std(firingRatesAverage(:,:,6000:10000),0,[2 3]));
z_score = (firingRatesAverage(:,:,:) - mu) ./ sigma;
mean_z_score = squeeze(mean(z_score(:,:,10000:14000), 3));
plot2 =subplot(1,2,2)
h2 = imagesc(mean_z_score);
caxis([-4  4]);
colormap(plot2,redblue);
colorbar
 
end