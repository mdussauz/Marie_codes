%plot average population PSTH 

firingRatesAverage = nanmean(smoothPSTH,5);
i = 1; 
x = -10:1.0000e-03:10-1.0000e-03;

figure()
for whatsmell = 1:5 
    for whatconc = 1:4
    subplot(5,4,i)
    PopAverage = squeeze(mean(firingRatesAverage(:,whatsmell,whatconc,:),1));
    p = plot(x,PopAverage);
    
    xlim([-9.8 9.8])
    ylim([0 14])
    p.LineWidth = 2;
    p.Color = [0 0.4470 0.7410];
    i = i+1;
    end
end

