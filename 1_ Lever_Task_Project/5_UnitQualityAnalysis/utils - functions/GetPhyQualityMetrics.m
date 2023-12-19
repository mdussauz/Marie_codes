%function [] = GetPhyQualityMetrics (SessionName)
%input:
% SessionName - eg =  'S12/2023-08-04_14-29-23'
% sessionHalfTime = Timepoint diving session in two to compare waveform
SessionName =  'S12/2023-08-04_14-29-23';
%% paths
if strcmp(computer, 'MACI64')
    ephyspath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Ephys/';
else
    ephyspath = '/mnt/data/Processed/Ephys/';
end

myKsDir = fullfile(ephyspath,SessionName);

%% Load data from kilosort/phy
if exist(fullfile(myKsDir, 'cluster_info.tsv')) 
   filename = fullfile(myKsDir, 'cluster_info.tsv');
end 

if ~isempty(filename)
    fid = fopen(filename);
    C = textscan(fid, '%s%s%s%s%s%s%s%s%s%s%s');
    fclose(fid);
end