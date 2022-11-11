% Replay_sniffing_20221108

SessionPath = 'O3/O3_20211005_r0_processed.mat';


if strcmp(computer, 'MACI64')
    datapath = '/Users/Priyanka/Desktop/LABWORK_II/Data/Smellocator/Processed/Behavior/';
elseif strcmp(computer,'PCWIN64')
    datapath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
else
    datapath = '/mnt/data/Processed/Behavior/';
end
MySession = fullfile(datapath,SessionPath);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

[TracesOut, ColNames, TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate] = LoadProcessedDataSession(MySession);

%% get sniffing
%% sniff trace
%sniffs = TracesOut(:,3);
fband = [0.5 30];
Np    = 4; % filter order
[b,a] = butter(Np,fband/(SampleRate/2)); % band pass Butterworth filter coefficients
TracesOut(:,3) = filtfilt(b,a,TracesOut(:,3)); %apply the filter to x(t)

% find peaks
[pks,pkloc] = findpeaks(TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);
% find valleys
[vls,vlloc] = findpeaks(-TracesOut(:,3),'MinPeakProminence',0.2,'MinPeakDistance',SampleRate*0.05);

% sanity checks
for i = 1:numel(pks)-1 
    % is there a unique valley
    thispkvalley = intersect(find(vlloc>pkloc(i)),find(vlloc<pkloc(i+1)));
    if numel(thispkvalley==1)
        Loc(i,:) = [pkloc(i) vlloc(thispkvalley)]; %sniff times
    else
        keyboard;
    end
end

%%
% figure; 
% subplot(1,3,1); 
% histogram(diff(Loc(:,1))); 
% subplot(1,3,2); 
% histogram(diff(Loc(:,2))); 
% subplot(1,3,3); 
% histogram((Loc(:,2)-Loc(:,1))); 

%% Trial Aligned Sniff times?
[AlignedSniffs] =  TrialAlignedSniffTimes(Loc,TracesOut(:,5),SampleRate);

%% Get the replay traces and spikes
MyUnits = 1; 
SingleUnits(1).spikes = Loc(:,1);
SingleUnits(1).trialtags = 0;
[MyTraces,timestamps,PSTH,Raster] = ...
ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
                        'plotfigures', 1, 'plotephys', 1, ...
                        'PlotOpenLoop', 1, ...
                            'whichunits', MyUnits);