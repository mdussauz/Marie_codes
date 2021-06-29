addpath(genpath('/opt/dPCA-master'))


NOdors = 20; %  for id/conc exp
NId = 5; % nber of odor id 
NConc = 4; % nber of od conc
NRep = 5; % nber of repeats 
Nneurons = size(smoothPSTH5D,1);

% 3 arrays needed to run code
NtrialNum = zeros(Nneurons, NId, NConc); %initialize 
NtrialNum(:,:,:) = 5;
trialNum = NtrialNum;
firingRates =  smoothPSTH5D;
firingRatesAverage = nanmean(firingRates,5);

N =  size(PSTH5D,1); % 97   % number of neurons
T = 20000;     % number of time points
S = 5;%7;       % number of odor
C = 4;          % number of concentration
E = 5;%20;     % maximal number of trial repetitions

time = (1:T) / 1000;
timeEvents = [10, 14];

% setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)


combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'identity', 'Concentration', 'Condition-independent', 'Id/Conc Interaction'};
margColours = [250 140 0; 0 115 255; 0 0 0; 150 150 150]/256; % for Florin

%time_window = [:]; % initial setting
%time_window = 1:20000; % to look at explained variance in a specific epoch of the trial
time_window = 10000:14000; % to look at explained variance in a specific epoch of the trial


%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage(:,:,:,time_window), W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);

%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage(:,:,:,time_window), W, V, ...
    'combinedParams', combinedParams);

% dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
%     'explainedVar', explVar, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg,                 ...
%     'time', time,                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3, ...
%     'legendSubplot',16);
%WITHOUT LEGEND
dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3);

