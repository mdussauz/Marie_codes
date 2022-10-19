function [] = PlotRespondTrend(smoothPSTH)

% written by MD % Plot trend of average response across concentrations
% input PSTH 4D from conc experiments
% Smoothened PSTH with the following dimensions:
% Neurons x Stimuli x time x Repeats 

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

%% sort units between excited and inhibited population
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
        odormean = squeeze(nanmean(firingRatesAverage(:,whatsmell,whatconc,10000:14000),[1 4]));
        normalized_mean = (odormean - airmean) / airmean;
        all_normalized_mean(whatconc,whatsmell) = normalized_mean;
         
        %FR Average for 1 stimulus for exc and inh populations 
        FRAverageExc = squeeze(excited_units(:,whatsmell,whatconc,:));
        FRAverageInh = squeeze(inhibited_units(:,whatsmell,whatconc,:));
        
        %Calculate FR Average during air and odor period for 1 stimulus for exc and inh populations 
        %Excitatory
        airmeanExc = squeeze(nanmean(excited_units(:,:,:,6000:10000),[1 2 3 4]));
        odormeanExc = squeeze(nanmean(excited_units(:,whatsmell,whatconc,10000:14000),[1 4]));
        %Inhibitory 
        airmeanInh = squeeze(nanmean(inhibited_units(:,:,:,6000:10000),[1 2 3 4]));
        odormeanInh = squeeze(nanmean(inhibited_units(:,whatsmell,whatconc,:),[1 4]));
        
        %Normalize FR averages 
        normalized_mean_Exc = (odormeanExc - airmeanExc) / airmeanExc;
        normalized_mean_Inh = (odormeanInh - airmeanInh) / airmeanInh;
        
        all_normalized_mean_Exc (whatconc,whatsmell) = normalized_mean_Exc;
        all_normalized_mean_Inh (whatconc,whatsmell) = normalized_mean_Inh;
    end
    
    X = [-4 -3 -2 -1]; %log concentration 
    Y_Exc = all_normalized_mean_Exc(:,whatsmell); %for each individual odor
    Y2_Exc = mean(all_normalized_mean_Exc,2); %mean across all odors
    Y_Inh = all_normalized_mean_Inh(:,whatsmell); %for each individual odor
    Y2_Inh = mean(all_normalized_mean_Inh,2); %mean across all odors
    Y_all = all_normalized_mean(:,whatsmell); %for each individual odor
    Y2_all = mean(all_normalized_mean,2); %mean across all odors;
    
    figure(1) %excitatory population plot
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
    title( ['All Odors'])
    
    figure(2) %inhibitory population plot
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
    title( ['All Odors'])
    
    figure (3) % entire population
    subplot(3,2,whatsmell)
    plot (X,Y_all, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(all_normalized_mean(:,whatsmell))-0.05 max(all_normalized_mean(:,whatsmell)+0.05)])
    xticks(X)
    yline(0,'k--');
    title( ['Odor',  num2str(whatsmell)])
    
    subplot(3,2,6)
    plot (X,Y2_all, '-b', 'LineWidth', 1.5)
    xlim([-5 0])
    ylim([min(Y2_all)-0.05 max(Y2_all)+0.05])
    xticks(X)
    yline(0,'k--');
    title( ['All Odors'])
    
end

end
