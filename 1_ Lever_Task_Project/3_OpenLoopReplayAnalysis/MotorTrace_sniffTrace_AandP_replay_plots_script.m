%Check motor and sniffing during active and passive replay
%Sanity check of replays

%% paths 
%classic open loop 
SessionName = 'O3/O3_20211005_r0_processed.mat';
% SessionName = 'O8/O8_20220702_r0_processed.mat';
% SessionName = 'O9/O9_20220630_r0_processed.mat';

%SessionName = 'S1/S1_20230314_r0_processed.mat';
%SessionName = 'S3/S3_20230321_r0_processed.mat';
%SessionName = 'S6/S6_20230727_r0_processed.mat';
%SessionName = 'S7/S7_20230707_r0_processed.mat';
%SessionName = 'S11/S11_20230812_r0_processed.mat';
%SessionName = 'S12/S12_20230727_r0_processed.mat';

%free lever open loop 
%SessionName = 'S1/S1_20230403_r0_processed.mat'; 
%SessionName = 'S6/S6_20230718_r0_processed.mat'; 
%SessionName = 'S7/S7_20230622_r0_processed.mat';
%SessionName = 'S11/S11_20230805_r0_processed.mat';

if strcmp(computer, 'MACI64')
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    ProcessedBehaviorPath = '/mnt/data/Processed/Behavior/';
end

MySession = fullfile(ProcessedBehaviorPath,SessionName);


%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');

%%
OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

%% Get the replay traces and spikes

NumActiveReplays = size(OpenLoop.ReplayTraces.Motor{1,1},2);
NumPassiveReplays = size(PassiveReplayTraces.Motor,2); 

figure(1) % Closed loop Template Motor and Sniffs Traces
subplot(2,1,1)
plot(OpenLoop.TemplateTraces.Motor{1,1})
title('Closed loop Template Motor Trace')
subplot(2,1,2)
plot(OpenLoop.TemplateTraces.Sniffs{1,1})
title('Closed loop Template Sniffs Trace')

figure(2) % Active Replay Motor Traces
subplot(NumActiveReplays+1,1,1)
plot(OpenLoop.ReplayTraces.Motor{1, 1}) ; %overlay of Active Replay Motor Traces
title('Overlay of Active Replay Motor Trace')

for i = 2:NumActiveReplays+1
    subplot(NumActiveReplays+1,1,i)
    plot(OpenLoop.ReplayTraces.Motor{1, 1}(:,i-1))
    if i ==2
        title('Active Replay Motor Trace for each repeat')
    end
end

figure(3) % Active Replay Sniff Traces
for j = 1:NumActiveReplays
    subplot(NumActiveReplays,1,j)
    plot(OpenLoop.ReplayTraces.Sniffs{1, 1}(:,j))
    if j ==1
        title('Active Replay Sniffs Trace for each repeat')
    end
end

figure(4) % Passive Replay Motor Traces
for i= 1:NumPassiveReplays
    subplot(NumPassiveReplays+1,1,1)
    plot(PassiveReplayTraces.Motor{1,i}); %overlay of Passive Replay Motor Traces
    if i ==1
        title('Overlay of Passive Replay Motor Traces')
    end
    hold on
    subplot(NumPassiveReplays+1,1,i+1)
    plot(PassiveReplayTraces.Motor{1,i});
        if i ==1
        title('Passive Replay Motor Trace for each repeat')
    end
end
hold off

figure(5) % Passive Replay Sniffs Traces
for j = 1:NumPassiveReplays
    subplot(NumPassiveReplays,1,j)
    plot(PassiveReplayTraces.Sniffs{1,j});
    if j ==1
    title('Passive Replay Sniffs Trace for each repeat')
    end
end


% for j = 2:NumActiveReplays+1
%     subplot(NumActiveReplays+1,1,j)
%     plot(OpenLoop.ReplayTraces.Motor{1, 1}(:,j-1))
% end

figure(6) %Active Replay Lever Trace
subplot(NumActiveReplays+1,1,1)
plot(OpenLoop.TemplateTraces.Lever{1,1})
title('Closed loop Lever Trace')
for i = 1:NumActiveReplays
    subplot(NumActiveReplays+1,1,i+1)
    plot(OpenLoop.ReplayTraces.Lever{1, 1}(:,i))
    if i ==1
        title('Active Replay Lever Trace for each repeat')
    end

y_cutoff = length(OpenLoop.ReplayTraces.Lever{1, 1}(:,1));
TemplateLeverTrace = OpenLoop.TemplateTraces.Lever{1,1}(1:y_cutoff);
%There are NaN in the lever trace so fill them
TemplateTrace = fillmissing(TemplateLeverTrace,'nearest');
distanceTemplate = sum(abs(diff(TemplateTrace)));
corr_activetrace(i) = corr(TemplateTrace,OpenLoop.ReplayTraces.Lever{1, 1}(:,i));
var_activetrace(i) = var(OpenLoop.ReplayTraces.Lever{1, 1}(:,i));
residuals_activetrace(i) = mean((TemplateTrace - OpenLoop.ReplayTraces.Lever{1, 1}(:,i)).^2)';
distance_ActiveReplays(i) = sum(abs(diff(OpenLoop.ReplayTraces.Lever{1, 1}(:,i))));

end

figure(7) %Passive Replay Lever Trace
subplot(NumPassiveReplays+1,1,1)
plot(OpenLoop.TemplateTraces.Lever{1,1})
title('Closed loop Lever Trace')
for i = 1:NumPassiveReplays
    subplot(NumPassiveReplays+1,1,i+1)
    plot(PassiveReplayTraces.Lever{1, i})
    if i ==1
        title('Active Replay Lever Trace for each repeat')
    end

    % in case templates and replay traces have different length:
    if length(PassiveReplayTraces.Lever{1, i}) < length(OpenLoop.TemplateTraces.Lever{1,1})
        y_cutoff = length(PassiveReplayTraces.Lever{1, i});
        TemplateLeverTrace = OpenLoop.TemplateTraces.Lever{1,1}(1:y_cutoff);
        ThisPassiveReplayLeverTrace = PassiveReplayTraces.Lever{1, i};
    elseif length(PassiveReplayTraces.Lever{1, i}) > length(OpenLoop.TemplateTraces.Lever{1,1})
        y_cutoff = length(PassiveReplayTraces.Lever{1, i});
        TemplateLeverTrace = OpenLoop.TemplateTraces.Lever{1,1};
        ThisPassiveReplayLeverTrace = PassiveReplayTraces.Lever{1, i}(1:y_cutoff);
    else
        TemplateLeverTrace = OpenLoop.TemplateTraces.Lever{1,1}(:);
        ThisPassiveReplayLeverTrace = PassiveReplayTraces.Lever{1, i};
    end
%If there are NaN in the lever trace, fill them
TemplateTrace = fillmissing(TemplateLeverTrace,'nearest');
distanceTemplate = sum(abs(diff(TemplateTrace)));
coor_passivetrace(i) = corr(TemplateTrace,ThisPassiveReplayLeverTrace);
var_passivetrace(i) = var(ThisPassiveReplayLeverTrace);
residuals_passivetrace(i) = mean((TemplateTrace - ThisPassiveReplayLeverTrace).^2)';
distance_PassiveReplays(i) = sum(abs(diff(ThisPassiveReplayLeverTrace)));

end



 %% Plot distance travelled with lever
distance = [distanceTemplate distance_ActiveReplays];

% rescaling distance travelled between 0 and 1 
minVal = min(distance);
maxVal = max(distance);
norm_data = (distance - minVal) / ( maxVal - minVal );
your_original_data = minVal + norm_data.*(maxVal - minVal); % just for sanity check

plot(norm_data, 'k-o')
xticks(2:2:10)
yticks(0:0.5:1)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)



