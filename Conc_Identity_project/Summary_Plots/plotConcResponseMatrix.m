function []= plotConcResponseMatrix(smoothPSTH)
%written by MD

firingRatesAverage = nanmean(smoothPSTH,5);
firstbin = -10;
lastbin = 10;

%% Reponse Matrix for each concentration and each odor
for whatsmell = 1:5
    figure(whatsmell)
    %figure()
    FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
    for whatconc = 1:4
        subplot(1,4,whatconc)
        X = FRAverageConc(:,whatconc,:);
        minvalue = min(X,[],'all');
        maxvalue = max(X,[],'all');
        allmin(whatsmell,whatconc) = minvalue;
        allmax(whatsmell,whatconc) = maxvalue;
        h = imagesc(squeeze(X));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10]);
        caxis([0  70]);
        %
        %        xline(0,'w--');
        %        xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear X; clear Y;
    end
    colormap(parula);
    title( ['Odor',  num2str(whatsmell)])
    
end


%% Z SCORE response matrix for each concentration and each odor
%
%for myNeuron = 1:Nneurons %for each neuron
for whatsmell = 1:5
    %figure()
    figure(whatsmell+5)
    for whatconc = 1:4
        FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
        subplot(1,4,whatconc)
        mu = squeeze(mean(FRAverageConc(:,:,6000:10000),[2 3]));
        sigma = squeeze(std(FRAverageConc(:,:,6000:10000),0,[2 3]));
        z_score = (FRAverageConc(:,whatconc,:) - mu) ./ sigma;
        h = imagesc(squeeze(z_score));
        set(h, 'XData', [firstbin, lastbin]);
        xlim([-10, 10])
        caxis([-4  4]);
        %         xline(0,'w--');
        %         xline(4,'w--');
        xline(0,'k--');
        xline(4,'k--');
        clear mu; clear sigma; clear z_score;
    end
    colormap(redblue);
    title( ['Odor',  num2str(whatsmell)])
    
end


%% Population mean reponse plots - trend of response
for whatsmell = 1:5
    %figure()
    for whatconc = 1:4
        FRAverageConc = squeeze(firingRatesAverage(:,whatsmell,:,:));
        airmean = squeeze(nanmean(FRAverageConc(:,:,6000:10000),[1 2 3]));
        odormean = squeeze(mean(FRAverageConc(:,whatconc,10000:14000),[1 3]));
        stdeviation = squeeze(std(FRAverageConc(:,whatconc,10000:14000),0,[1 3]));
        normalized_mean = (odormean - airmean) / airmean;
        all_normalized_mean (whatconc,whatsmell) = normalized_mean;
    end
    
    X = [-4 -3 -2 -1];
    Y = all_normalized_mean(:,whatsmell);
    Y2 = mean(all_normalized_mean,2);
    
    figure(11)
    subplot(3,2,whatsmell)
    plot (X,Y, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(all_normalized_mean(:,whatsmell))-0.05 max(all_normalized_mean(:,whatsmell)+0.05)])
    xticks(X)
    yline(0,'k--');
    title( ['Odor',  num2str(whatsmell)])
    
    subplot(3,2,6)
    plot (X,Y2, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(Y2)-0.05 max(Y2)+0.05])
    xticks(X)
    yline(0,'k--');
    
    
end

%% Mean response for each repeat (column) - variability across repeats

for whatsmell = 1:5
    figure(whatsmell+11)
    for whatconc = 1:4
        subplot(1,4,whatconc)
        for rep = 1:5
            FRAverageRep = squeeze(smoothPSTH(:,whatsmell,whatconc,:,rep));
            odormeanRep = squeeze(mean(FRAverageRep(:,10000:14000),2));
            all_rep_mean (:,rep) = odormeanRep;
        end
        h = imagesc(all_rep_mean);
        caxis([0  50]); %80
        colormap(redblue);
        clear all_rep_mean
    end
    title( ['Odor',  num2str(whatsmell)])
end


end