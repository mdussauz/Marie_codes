% Sort neurons into responsive excitatory and inhibitory neurons

firingRatesAverage = nanmean(smoothPSTH,5);
Nneurons = size(smoothPSTH, 1);

%% perform t test to find responsive neurons
for n = 1:Nneurons %for each neuron
    air = squeeze(firingRatesAverage(n,:,:,6000:10000));
    odor = squeeze(firingRatesAverage(n,:,:,10000:14000));
    t_test(n,:,:) = ttest (odor, air, 'alpha', 0.05,'Dim', 3);
    
end

%% sort responsive units between excited and inhibited
for whatsmell = 1:5
    for whatconc = 1:4
        FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
        mu = squeeze(mean(FRAverageConc(:,:,1:10000),[2 3]));
        sigma = squeeze(std(FRAverageConc(:,:,1:10000),0,[2 3]));
        z_score = (FRAverageConc(:,whatconc,:) - mu) ./ sigma;
        z_score_matrix(:,whatsmell,whatconc) = mean(z_score(:,10000:14000),2);
        clear mu; clear sigma; clear z_score;
    end
end

responsive = t_test.*z_score_matrix;

excited_units = zeros(size(firingRatesAverage));
inhibited_units = zeros(size(firingRatesAverage));

for n = 1:Nneurons
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
    %figure()
    for whatconc = 1:4
        FRAverageConcExc = squeeze(excited_units(:,whatsmell,:,:));
        FRAverageConcInh = squeeze(inhibited_units(:,whatsmell,:,:));
        airmeanExc = squeeze(nanmean(FRAverageConcExc(:,:,6000:10000),[1 2 3]));
        odormeanExc = squeeze(nanmean(FRAverageConcExc(:,whatconc,10000:14000),[1 3]));
        airmeanInh = squeeze(nanmean(FRAverageConcInh(:,:,6000:10000),[1 2 3]));
        odormeanInh = squeeze(nanmean(FRAverageConcInh(:,whatconc,10000:14000),[1 3]));
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

%% plot population average response in time
i = 1;
x = -10:1.0000e-03:10-1.0000e-03;
figure(3)
for whatsmell = 1:5
    for whatconc = 1:4
        subplot(5,4,i)
        PopAverage = squeeze(nanmean(excited_units(:,whatsmell,whatconc,:),1));
        p = plot(x,PopAverage);
        
        xlim([-9.8 9.8])
        ylim([0 20])
        p.LineWidth = 2;
        p.Color = [0 0.4470 0.7410];
        i = i+1;
    end
end

j = 1;
figure(4)
for whatsmell = 1:5
    for whatconc = 1:4
        subplot(5,4,j)
        PopAverage = squeeze(nanmean(inhibited_units(:,whatsmell,whatconc,:),1));
        p = plot(x,PopAverage);
        
        xlim([-9.8 9.8])
        ylim([0 7])
        p.LineWidth = 2;
        p.Color = [0 0.4470 0.7410];
        j = j+1;
    end
end