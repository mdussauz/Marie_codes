function [W, explVar,whichMarg] = RundPCA(PSTH5D)
addpath(genpath('/opt/dPCA-master'))

% dPCA 
%The data should be joined in three arrays of the following sizes :
%
% trialNum: N x S x D
% firingRates: N x S x D x T x maxTrialNum
% firingRatesAverage: N x S x D x T
%
% N is the number of neurons
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
% D is the number of decisions (D=2)
% T is the number of time-points (note that all the trials should have the
% same length in time!)
%
% trialNum -- number of trials for each neuron in each S,D condition (is
% usually different for different conditions and different sessions)
%
% firingRates -- all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.
%
% firingRatesAverage -- average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates,5)
% If it's filled up with zeros (as is convenient if it's stored on hard 
% drive as a sparse matrix), then 
%    firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRates,5)./trialNum)

decoding = 0; % optional 


Nneurons = size(PSTH5D,1);
%NOdors = size(PSTH5D,2)*size(PSTH5D,3); %  for id/conc exp % unused
NId = size(PSTH5D,2); % nber of odor id 
NConc = size(PSTH5D,3); % nber of od conc
NRep = size(PSTH5D,5); % nber of repeats 

% 3 arrays needed to run code
NtrialNum = zeros(Nneurons, NId, NConc); %initialize 
NtrialNum(:,:,:) = NRep; %used to be 5
trialNum = NtrialNum;
firingRates =  PSTH5D;
firingRatesAverage = nanmean(firingRates,5);

N =  size(PSTH5D,1); % 97   % number of neurons
T = size(PSTH5D,4);     % number of time points
S = size(PSTH5D,2);%7;       % number of odor
C = size(PSTH5D,3);          % number of concentration
E = size(PSTH5D,5);%20;     % maximal number of trial repetitions

time = (1:T) / 1000;
timeEvents = [6, 8]; %time of odor

% setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)

%we have the following parameters:
%    1 - odor 
%    2 - concentration
%    3 - time
% There are three pairwise interactions:
%    [1 3] - odor/time interaction
%    [2 3] - concentration/time interaction
%    [1 2] - odor/concentration
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'identity', 'Concentration', 'Condition-independent', 'Id/Conc Interaction'};
%margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
margColours = [250 140 0; 0 115 255; 0 0 0; 150 150 150]/256; % for Florin


% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 
%    combinedParams = {{1, [1 2]}, {2}};

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
%timeEvents = time(round(length(time)/2));
%timeEvents = time(round(length(time)/2));

% check consistency between trialNum and firingRates
% for n = 1:size(firingRates,1)
%     for s = 1:size(firingRates,2)
%         for d = 1:size(firingRates,3)
%             assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
%         end
%     end
% end

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
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

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
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


%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3);

%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: decoding

if decoding == 1
    decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};
    
    accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 5, ...        % increase to 100
        'filename', 'tmp_classification_accuracy.mat');
    
    dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
    
    accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
        'lambda', optimalLambda, ...
        'combinedParams', combinedParams, ...
        'decodingClasses', decodingClasses, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 5, ...        % increase to 100
        'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
        'filename', 'tmp_classification_accuracy.mat');
    
    dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
    
    componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);
    
    dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
        'explainedVar', explVar, ...
        'marginalizationNames', margNames, ...
        'marginalizationColours', margColours, ...
        'whichMarg', whichMarg,                 ...
        'time', time,                        ...
        'timeEvents', timeEvents,               ...
        'timeMarginalization', 3,           ...
        'legendSubplot', 16,                ...
        'componentsSignif', componentsSignif);
end

end
