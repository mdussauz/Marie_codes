%Check motor and sniffing during active and passive replay

%classic open loop 
MySession = '/mnt/grid-hs/pgupta/Behavior/S1/S1_20230314_r0.mat';
% PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/S3/S3_20230321_r0.mat',1);
% PreprocessSmellocatorData('mnt/grid-hs/pgupta/Behavior/S6/S6_20230727_r0.mat',1);
% PreprocessSmellocatorData('mnt/grid-hs/pgupta/Behavior/S7/S7_20230707_r0.mat',1);
% PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/S11/S11_20230805_r0.mat',1);
% PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/S12/S12_20230727_r0.mat',1);

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% Get the replay traces and spikes