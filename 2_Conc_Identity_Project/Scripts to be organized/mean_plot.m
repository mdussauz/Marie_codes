for whatconc = 1:4
        airmean = squeeze(nanmean(FRAverageConc(:,:,6000:10000),[1 2 3]));
        odormean = squeeze(mean(FRAverageConc(:,whatconc,10000:14000),[1 3]));
        stdeviation = squeeze(std(FRAverageConc(:,whatconc,10000:14000),0,[1 3]));
        normalized_mean = (odormean - airmean) / airmean;
        all_normalized_mean (whatconc) = normalized_mean;
end

X = [-4 -3 -2 -1];
Y = all_normalized_mean;

figure(4)
plot (X,Y, '-b', 'LineWidth', 1.5)
xlim([-5 0])
ylim([0.02 0.1])
xticks(X)
yline(0,'k--');