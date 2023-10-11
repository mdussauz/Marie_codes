close all
clearvars

% Change path to the folder you saved the data I shared with you
addpath     'C:\Users\deher\Dropbox\Data_Marie\Responsiveness_test'
load        'C:\Users\deher\Dropbox\Data_Marie\Responsiveness_test\test_dataset.mat' 

% This code classifies unit/channel/bouton/neuron responsiveness by
% evaluating if certain ranges of the trace (frame windows) are over
% (enhanced) or under (suppressed) the defined percentile threshold of the
% baseline distribution values. As the method assumes that the baseline has
% a normal distribution, you can check the baseline distributions for each
% channel if you leave 'baseline_plots = 1'.

%% Basic variables needed to run the code

reference                   = [56 62; 63 69; 70 76; 77 83; 84 90; 91 97];  % I am using test windows with 7 frames starting at the cue delivery period. It is important to use the same window size for each window you want to evaluate.
reference_baseline          = [28 34; 35 41; 42 48; 49 55];                % You have to divide your baseline in frame windows equivalent to the ones you want to evaluate.
reference_threshold         = [1 99; 5 95; 10 90];                         % Here you can feed possible thresholds you want to try. To get an idea of which threshold to pick, use 'Responsive_bouton_threshold_test_20221021.m'                                   
threshold                   = 1;                                           % Pick the desired threshold by defining the used row of reference threshold.
baseline_plots              = 1;                                           % 1: ON / 0: OFF option for baseline distribution plots for each channel.                                                                   
perc_plots                  = 1;                                           % 1: ON / 0: OFF option for percentage analysis plots.                    

%% Threshold test - use this to check which threshold to pick and use in 'reference_threshold'

[f1, resptest_perctest]... 
    = Responsive_bouton_threshold_test_20221021(...
        dataset, reference, reference_baseline);

f1.Position = [550 300 800 500];


%% Responsiveness evaluation - use this to get the results

[f2, f3, resptest_response_class]... 
    = Responsive_bouton_ID_20221021(...
        dataset, reference, reference_baseline, reference_threshold, ...
        threshold, baseline_plots, perc_plots);

maxfig(f2,1);
f3.Position = [600 350 800 500];

% The results for each unit/channel/bouton/neuron responsiveness are in resptest_response_class: 
% 0: Unresponsive
% 1: Enhanced
% 2: Suppressed
% 3: Complex