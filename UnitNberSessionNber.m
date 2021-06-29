%% USER: what brain area

brain_area = 'AON'; % input is 'AON' or 'APC'

%% add the relevant repositories to path
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
addpath(genpath('/opt/afterphy'))
addpath(genpath('/opt/spikes'))
addpath(genpath('/opt/npy-matlab'));
addpath(genpath('/opt/Marie_codes'));
addpath(genpath('/opt/Marie_codes/dPCA'))

%% Filepaths
%folder with stimulus files
addpath(genpath('/mnt/data/PhotonCerber_Stimuli'))
%folders with sessions coming from AON mice
addpath (genpath('/mnt/interim/N01'))
addpath (genpath('/mnt/interim/N02'))
addpath (genpath('/mnt/interim/D1'))
addpath (genpath('/mnt/interim/D2'))
addpath (genpath('/mnt/interim/D3'))
%folders with sessions coming from APC mice
addpath (genpath('/mnt/interim/PJ1'))
addpath (genpath('/mnt/interim/PJ2'))
addpath (genpath('/mnt/interim/APC1'))
addpath (genpath('/mnt/interim/APC2'))

%% Recording sessions

if brain_area == 'AON'
    myKsDir = {'/mnt/interim/N01/2021-01-25_09-35-25','/mnt/interim/N01/2021-01-27_10-46-51',...
        '/mnt/interim/N02/2021-01-26_11-47-01','/mnt/interim/N02/2021-01-28_10-57-31',...
        '/mnt/interim/N02/2021-01-29_14-06-12',...
        '/mnt/interim/D1/2021-04-27_09-30-31', '/mnt/interim/D1/2021-05-05_11-42-47',...
        '/mnt/interim/D1/2021-05-06_13-53-52', '/mnt/interim/D1/2021-05-07_09-13-39',...
        '/mnt/interim/D3/2021-04-26_17-24-41', '/mnt/interim/D3/2021-04-28_15-59-15',...
        '/mnt/interim/D3/2021-04-29_16-43-06','/mnt/interim/D3/2021-04-30_16-23-08',...
        '/mnt/interim/D3/2021-05-20_13-54-53','/mnt/interim/D3/2021-05-21_13-04-28',...
        '/mnt/interim/D3/2021-05-24_14-33-50','/mnt/interim/D3/2021-05-25_14-39-17',...
        '/mnt/interim/D3/2021-05-26_13-18-29','/mnt/interim/D3/2021-05-27_14-11-45',...
        '/mnt/interim/D2/2021-05-31_13-19-08', '/mnt/interim/D2/2021-06-01_12-29-35',...
        '/mnt/interim/D2/2021-06-02_13-28-26', '/mnt/interim/D2/2021-06-03_15-09-14',...
        '/mnt/interim/D2/2021-06-09_12-39-32','/mnt/interim/D2/2021-06-10_12-56-18'...
        };
elseif brain_area == 'APC'
    myKsDir = {'/mnt/interim/PJ1/2021-01-25_15-13-27','/mnt/interim/PJ1/2021-01-28_14-21-40',...
        '/mnt/interim/PJ2/2021-02-14_15-06-20','/mnt/interim/PJ2/2021-02-15_14-51-10',...
        '/mnt/interim/APC1/2021-04-26_11-12-59', '/mnt/interim/APC1/2021-04-27_15-37-37',...
        '/mnt/interim/APC1/2021-04-28_09-23-05','/mnt/interim/APC1/2021-04-29_10-12-26',...
        '/mnt/interim/APC1/2021-05-03_14-58-49', '/mnt/interim/APC1/2021-05-04_14-16-04',...
        '/mnt/interim/APC2/2021-04-27_12-35-49', '/mnt/interim/APC2/2021-04-28_12-09-32',...
        '/mnt/interim/APC2/2021-04-29_13-19-18', '/mnt/interim/APC2/2021-04-30_13-26-00',...
        '/mnt/interim/APC2/2021-05-03_11-12-06', '/mnt/interim/APC2/2021-05-04_11-24-28'...
        };
end

%% Stimulus File

if brain_area == 'AON'
    stimfilename = {'210125_9_35.txt', '210127_10_47.txt',... %N01
        '210126_11_47.txt','210128_10_58.txt', '210129_14_06.txt',...%N02
        '210427_9_30.txt','21055_11_42.txt','21056_13_53.txt','21057_9_13.txt',... %D1
        '210426_17_24.txt','210428_15_58.txt','210429_16_42.txt','210430_16_22.txt',...%D3
        '210520_13_55.txt', '210521_13_04.txt', '210524_14_33.txt',...
        '210525_14_39.txt','210526_13_18.txt', '210527_14_11.txt',...
        '210531_13_19.txt', '21061_12_29.txt', '21062_13_33.txt', '21063_15_09.txt', ...%D2 %'21062_13_28.txt' is incorrect
        '21069_12_39.txt','210610_12_56.txt'...
        };
elseif brain_area == 'APC'
    stimfilename = {'210125_15_13.txt', '210128_14_21.txt',... %PJ1
        '210214_15_06.txt', '210215_14_52.txt',... %PJ2
        '210426_11_12.txt', '210427_15_37.txt', '210428_9_22.txt',...%APC1
        '210429_10_11.txt','210430_9_23.txt','21053_14_58.txt','21054_14_15.txt',...
        '210427_12_35.txt', '210428_12_09.txt','210429_13_18.txt', ... %APC2
        '210430_13_25.txt','21053_11_11.txt','21054_11_23.txt',...
        };
end



%% Get the structures  

for i = 1: length(myKsDir) % loop through each recording session

    
    %% Load data from kilosort/phy
    sp = loadKSdir(myKsDir{i}); % get all units
    
    %% Number of good units 
    
    goodUNber = length (find(sp.cgs==2)); % count how many good units
    disp(['found ',num2str(goodUNber),' good units']); % display for user
    
    UnitsPerSession(i) =  goodUNber;
    
end 

NberOfSessions = length(UnitsPerSession);