%Check motor and sniffing during active and passive replay
%Sanity check of replays

%classic open loop 
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S1/S1_20230314_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S3/S3_20230321_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S6/S6_20230727_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S7/S7_20230707_r0_processed.mat';
MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S11/S11_20230812_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S12/S12_20230727_r0_processed.mat';

%free lever open loop 
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S1/S1_20230403_r0_processed.mat'; 
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S6/S6_20230718_r0_processed.mat'; 
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S7/S7_20230622_r0_processed.mat';
%MySession = '/mnt/grid-hs/mdussauz/Smellocator/Processed/Behavior/S11/S11_20230805_r0_processed.mat';


%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% Get the replay traces and spikes

NumActiveReplays = size(OpenLoop.ReplayTraces.Motor{1,1},2);
NumPassiveReplays = size(PassiveReplayTraces.Motor,2); 

figure(1) % Closed loop Template Motor and Sniffs Traces
title('Closed loop Template Motor and Sniffs Trace')
subplot(2,1,1)
plot(OpenLoop.TemplateTraces.Motor{1,1})
subplot(2,1,2)
plot(OpenLoop.TemplateTraces.Sniffs{1,1})

figure(2) % Active Replay Motor Traces
title('Active Replay Motor Trace')
subplot(NumActiveReplays+1,1,1)
plot(OpenLoop.ReplayTraces.Motor{1, 1}) ; %overlay of Active Replay Motor Traces

for i = 2:NumActiveReplays+1
    subplot(NumActiveReplays+1,1,i)
    plot(OpenLoop.ReplayTraces.Motor{1, 1}(:,i-1))
end

figure(3) % Active Replay Sniff Traces
title('Active Replay Sniffs Trace')
for j = 1:NumActiveReplays
    subplot(NumActiveReplays,1,j)
    plot(OpenLoop.ReplayTraces.Sniffs{1, 1}(:,j))
end

figure(4) % Passive Replay Motor Traces
title('Passive Replay Motor Trace')
for i= 1:NumPassiveReplays
    subplot(NumPassiveReplays+1,1,1)
    plot(PassiveReplayTraces.Motor{1,i}); %overlay of Passive Replay Motor Traces
    hold on
    subplot(NumPassiveReplays+1,1,i+1)
    plot(PassiveReplayTraces.Motor{1,i});
end
hold off

figure(5) % Passive Replay Sniffs Traces
title('Passive Replay Sniffs Trace')
for j = 1:NumPassiveReplays
    subplot(NumPassiveReplays,1,j)
    plot(PassiveReplayTraces.Sniffs{1,j});
end


% for j = 2:NumActiveReplays+1
%     subplot(NumActiveReplays+1,1,j)
%     plot(OpenLoop.ReplayTraces.Motor{1, 1}(:,j-1))
% end


