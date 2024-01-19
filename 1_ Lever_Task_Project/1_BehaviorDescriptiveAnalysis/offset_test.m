%S1 and S6 bad

%load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S1/S1_20230405_r0_processed.mat')
%load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S3/S3_20230405_r0_processed.mat')
%load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S6/S6_20230721_r0_processed.mat')
%load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S7/S7_20230616_r0_processed.mat')
load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S11/S11_20230807_r0_processed.mat')
%load('/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/S12/S12_20230806_r0_processed.mat')
%% not sorted

PlotLocationOffsetTracesMD(Traces,TrialInfo)
% whichcluster = 1;
% spiketimes = SingleUnits(whichcluster).spikes;
% PlotLocationOffsetPSTH_2018(whichcluster,Traces,TrialInfo,spiketimes)