function [handles] = LoadKilosortDefaultsMD(handles, Username)

KiloSortPath = '/opt/KiloSort/';
addpath(genpath(KiloSortPath)) % path to kilosort folder
addpath(genpath('/opt/npy-matlab/')) % path to npy-matlab scripts

% default settings
handles.FilePaths.Data(1) = {'/mnt/grid-hs/pgupta/EphysData'}; % Root storage
handles.FilePaths.Data(2) = {'PCX4'}; % local read/write folder
handles.FilePaths.Data(3) = {'/mnt/data/Sorted/'}; % local read/write folder
handles.ServerPath = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Ephys';
handles.YourConfigFile = fullfile(KiloSortPath,'StandardConfig_Albeanu.m');

% spike detection settings
handles.init_from_data = 0; % generate template spikes from data
handles.spike_det_settings.Data(1) = 2; % x, number of clusters - x times more than Nchan
handles.spike_det_settings.Data(2) = -4; % spike threshold in standard deviations (4)
handles.filter2binary.Value = 1; % default is to save bandpassed filtered data locally to binary file
handles.computeCAR.Value = 1; % default is to subtract the Common average reference
handles.IgnoreChannels.String = ''; % channels that shouldn't be included in computing CAR

% recording settings
handles.recording_settings.Data(1) = 64; % no. of spike channels saved
handles.InactiveChannels.String = ''; % channels that shouldn't be loaded

% reorder channels - for the new EIB
handles.ReorderChannels = [];

handles.auxchannels = []; % only used for binary files

% overwrite settings as per need
switch Username
    case {'K4'}
        handles.FilePaths.Data(1) = {'/mnt/data/Priyanka'}; % Root storage
        handles.FilePaths.Data(2) = {'K4'}; % Animal Name
        handles.IgnoreChannels.String = mat2str([22 26 41 42 43 45 46 49 50 53 54 55 56 59 60]);
    case {'PCX1','PCX3','PCX4','PCX5'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/pgupta/EphysData'};
        % handles.FilePaths.Data(1) = {'/mnt/data/EphysRaw'};
        handles.FilePaths.Data(2) = {Username};
        handles.recording_settings.Data(1) = 64;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB64');
        handles.ReorderChannels = EIB64;
    case {'APC1'}
        handles.FilePaths.Data(1) = {'/mnt/data/Sorted'};
        handles.FilePaths.Data(2) = {Username};
        handles.recording_settings.Data(1) = 64;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([61:64]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB64');
        handles.ReorderChannels = EIB64;
    case {'APC2'}
        handles.FilePaths.Data(1) = {'/mnt/data/Sorted'};
        handles.FilePaths.Data(2) = {Username};
        handles.recording_settings.Data(1) = 64;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([57:64]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB64');
        handles.ReorderChannels = EIB64;   
    case {'APC1binary'}
        handles.FilePaths.Data(1) = {'/mnt/data/Sorted'};
        handles.FilePaths.Data(2) = {'APC1'};
        handles.recording_settings.Data(1) = 72;
        handles.auxchannels = 0;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([61:72]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB64');
        handles.ReorderChannels = horzcat(EIB64, 65:72);
    case {'O5','O3'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/mdussauz/ephysdata/lever_task/BatchO'};
        %handles.FilePaths.Data(1) = {'/mnt/data/EphysRaw'};
        handles.FilePaths.Data(2) = {Username};
        handles.recording_settings.Data(1) = 64;
        handles.InactiveChannels.String = mat2str(9:32);
        handles.IgnoreChannels.String = mat2str([]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB64');
        handles.ReorderChannels = EIB64;
    case {'O2','O1'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/mdussauz/ephysdata/lever_task/BatchO'};
        %handles.FilePaths.Data(1) = {'/mnt/data/EphysRaw'};
        handles.FilePaths.Data(2) = {Username};
        handles.recording_settings.Data(1) = 64;
        handles.InactiveChannels.String = mat2str(1:32);
        handles.IgnoreChannels.String = mat2str([]);
        load(fullfile(KiloSortPath,'Priyanka','EIB_maps.mat'),'EIB32_new');
        handles.ReorderChannels = horzcat(EIB32_new,EIB32_new + 32);
    case {'MO1'}
        handles.FilePaths.Data(1) = {'/mnt/data/Priyanka'};
        handles.FilePaths.Data(2) = {'MO1'};
        handles.recording_settings.Data(1) = 32;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);
    case {'Batch_E_2_3'}
        %AON: E2, E3; APC: E6
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp'};
        handles.FilePaths.Data(2) = {'E2'};
        handles.recording_settings.Data(1) = 32;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);

    case {'Batch_E_6'}
        %AON: E2, E3; APC: E6
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/mdussauz/ephysdata/Conc_id_exp'};
        handles.FilePaths.Data(2) = {'E6'};
        handles.recording_settings.Data(1) = 40;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);

    case {'Batch_O6_9'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/mdussauz/ephysdata/lever_task/BatchO'};
        handles.FilePaths.Data(2) = {'O6'}; % change name of mouse here
        handles.recording_settings.Data(1) = 40;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);

    case {'Blom_Behavior_Mice'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/pgupta/EphysData/odor_location'};
        handles.FilePaths.Data(2) = {'O6'}; % change name of mouse here
        handles.recording_settings.Data(1) = 40;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);
    case {'Blom_Conc_Mice'}
        handles.FilePaths.Data(1) = {'/mnt/grid-hs/pgupta/EphysData/odor_location'};
        handles.FilePaths.Data(2) = {'E2'}; % change name of mouse here: E2 or E3
        handles.recording_settings.Data(1) = 32;
        handles.InactiveChannels.String = mat2str([]);
        handles.IgnoreChannels.String = mat2str([]);
    
end
end