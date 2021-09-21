function [] = PlotRespondTrend(smoothPSTH)

tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]
firingRatesAverage = nanmean(x,5);

%% calculate z score
mu = squeeze(mean(firingRatesAverage(:,:,:,6000:10000),[2 3 4]));
sigma = squeeze(std(firingRatesAverage(:,:,:,6000:10000),0,[2 3 4]));
z_score = (firingRatesAverage(:,:,:,:) - mu) ./ sigma;
z_score_matrix(:,:,:) = mean(z_score(:,:,:,10000:14000),4);

%% sort responsive units between excited and inhibited population
excited_units = zeros(size(firingRatesAverage));
inhibited_units = zeros(size(firingRatesAverage));

for n = 1:ncells
    for whatsmell = 1:5
        for whatconc = 1:4
            if z_score_matrix(n,whatsmell,whatconc) >0
                excited_units(n,whatsmell,whatconc,:) = firingRatesAverage(n,whatsmell,whatconc,:);
                inhibited_units(n,whatsmell,whatconc,:) = NaN;
            elseif z_score_matrix(n,whatsmell,whatconc) <0
                inhibited_units(n,whatsmell,whatconc,:) = firingRatesAverage(n,whatsmell,whatconc,:);
                excited_units(n,whatsmell,whatconc,:) = NaN;
            else
                excited_units(n,whatsmell,whatconc,:) = NaN;
                inhibited_units(n,whatsmell,whatconc,:) = NaN;
            end
        end
    end
end

%% plot mean population response trend
for whatsmell = 1:5
    
    for whatconc = 1:4
        %Normalized mean for the whole population
        airmean = squeeze(nanmean(firingRatesAverage(:,:,:,6000:10000),[1 2 3 4]));
        odormean = squeeze(nanmean(firingRatesAverage(:,whatsmell,whatconc,10000:14000),[1 3]));
        normalized_mean = (odormean - airmean) / airmean;
        
        %FR Average for 1 stimulus for exc and inh populations 
        FRAverageExc = squeeze(excited_units(:,whatsmell,whatconc,:));
        FRAverageInh = squeeze(inhibited_units(:,whatsmell,whatconc,:));
        
        %Calculate FR Average during air and odor period for 1 stimulus for exc and inh populations 
        %Excitatory
        airmeanExc = squeeze(nanmean(FRAverageConcExc(:,:,6000:10000),[1 2 3]));
        odormeanExc = squeeze(nanmean(FRAverageConcExc(:,whatconc,10000:14000),[1 3]));
        %Inhibitory 
        airmeanInh = squeeze(nanmean(FRAverageConcInh(:,:,6000:10000),[1 2 3]));
        odormeanInh = squeeze(nanmean(FRAverageConcInh(:,whatconc,10000:14000),[1 3]));
        
        %Normalize FR averages 
        normalized_mean_Exc = (odormeanExc - airmeanExc) / airmeanExc;
        normalized_mean_Inh = (odormeanInh - airmeanInh) / airmeanInh;
        
        all_normalized_mean_Exc (whatconc,whatsmell) = normalized_mean_Exc;
        all_normalized_mean_Inh (whatconc,whatsmell) = normalized_mean_Inh;
    end
    
    X = [-4 -3 -2 -1];
    Y_Exc = all_normalized_mean_Exc(:,whatsmell);
    Y2_Exc = mean(all_normalized_mean_Exc,2);
    Y_Inh = all_normalized_mean_Inh(:,whatsmell);
    Y2_Inh = mean(all_normalized_mean_Inh,2);
    
    figure(1)
    subplot(3,2,whatsmell)
    plot (X,Y_Exc, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(all_normalized_mean_Exc(:,whatsmell))-0.05 max(all_normalized_mean_Exc(:,whatsmell)+0.05)])
    xticks(X)
    yline(0,'k--');
    title( ['Odor',  num2str(whatsmell)])
    
    subplot(3,2,6)
    plot (X,Y2_Exc, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(Y2_Exc)-0.05 max(Y2_Exc)+0.05])
    xticks(X)
    yline(0,'k--');
    
    figure(2)
    subplot(3,2,whatsmell)
    plot (X,Y_Inh, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(all_normalized_mean_Inh(:,whatsmell))-0.05 max(all_normalized_mean_Inh(:,whatsmell)+0.05)])
    xticks(X)
    yline(0,'k--');
    title( ['Odor',  num2str(whatsmell)])
    
    subplot(3,2,6)
    plot (X,Y2_Inh, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(Y2_Inh)-0.05 max(Y2_Inh)+0.05])
    xticks(X)
    yline(0,'k--');
    
    
end

end
