function[]= plotSpikesinPhase(cluster, ind, trialInd, phaseStarts, resp)
% resp includes timebase, respiration_ds_filt, idx, peaks, plot(boolean)
spikes = cluster.spikes{ind, trialInd};
odorspikes = cluster.odorspikes{ind, trialInd};
y = zeros(size(spikes));
for ind = 1:length(spikes)
    y(ind) = getPolar(spikes(ind), phaseStarts);
end



figure();
subplot(4,3,[1 6]);
scatter(odorspikes, y)
xline(0)
xline(4)
ylabel("Phase");
set(gca, 'XLim', [-10, 10]);

subplot(4, 3, [7 9]);
plot(1:length(resp.Respiration_DS_Filt),resp.Respiration_DS_Filt);
hold on
plot(resp.Idx,-resp.MyPeaks,'or');
set(gca,'XLim',[findClosestInd(spikes(1), resp.newTimeBase), findClosestInd(spikes(length(spikes)), resp.newTimeBase)]);
hold off
ylabel("Respiration");


subplot(4, 3, [10 12]);
firstbin = -10;
lastbin = 10;
step = 0.1;
tbins = firstbin:step:lastbin;
myspikecount = histcounts(odorspikes, tbins).*(1/step);
smoothcount = smoothdata(myspikecount, 'gaussian', 10);
plot(smoothcount)
ylabel("Firing Rate");
xlabel("Time");

end



    
function[ind] = findClosestInd(time, timebase)
    [~, ind] = min(abs(timebase - time));
end




