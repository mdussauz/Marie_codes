function [] = PlotFiringRateDistribution(smoothPSTH)

%written by MD
%visualize distribution of mean FR during baseline and odor epoch 
%input PSTH 4D from conc experiments
% Smoothened PSTH with the following dimensions:
% Neurons x Stimuli x time x Repeats 

tot_time = 20000;
ncells = size (smoothPSTH,1);
max_rep = 5;
x = reshape(smoothPSTH,[ncells 5 4 tot_time max_rep]); %[nber of cells id conc timestamps rep]

%% plot FR distribution during baseline epoch
baselineFR = squeeze(mean(x(:,:,:,6000:10000,:),4));
baselineFR_vector = baselineFR(:); % reshaping into a vector
figure();
subplot(2,2,1); histogram(baselineFR_vector, 'BinMethod','integers','DisplayStyle','stairs')
title ('baseline FR distribution');
subplot(2,2,2); boxplot(baselineFR_vector,'PlotStyle','compact','Whisker',Inf,'orientation', 'horizontal')
%subplot(2,2,2); boxplot(baselineFR_vector,'PlotStyle','compact','Symbol','','orientation', 'horizontal') % to hide outliers

%% plot FR distribution during odor epoch
odorFR = squeeze(mean(x(:,:,:,10000:14000,:),4));
odorFR_vector = odorFR(:);% reshaping into a vector
subplot(2,2,3); histogram(odorFR_vector, 'BinMethod','integers','DisplayStyle','stairs')
title ('odor FR distribution');
subplot(2,2,4); boxplot(odorFR_vector,'PlotStyle','compact','Whisker',Inf,'orientation', 'horizontal') % to hide outliers
%subplot(2,2,4); boxplot(odorFR_vector,'PlotStyle','compact','Symbol','',
%'orientation', 'horizontal') % to hide outliers

end 