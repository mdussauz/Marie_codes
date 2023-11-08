addpath(genpath('/Users/Priyanka/Desktop/github_local/open-ephys'))

%% read and process event timestamps from open ephys file
[Events] = ParseOpenEphysEvents('all_channels.events');
% trial events are in Channel1, odor events are in Channel2

TrialTimeStamps = [Events.Channel4.On Events.Channel4.Off];
% delete first entry - just an empty trigger
TrialTimeStamps(1,:) = [];
OdorTimeStamps = [Events.Channel1.On Events.Channel1.Off];


% read PID data from open ephys file
[PhotoReceiver1, Timestamps, info] = load_open_ephys_data('100_ADC3.continuous');
[PhotoReceiver2, Timestamps, info] = load_open_ephys_data('100_ADC4.continuous');
[LED1, Timestamps, info] = load_open_ephys_data('100_ADC5.continuous');
[LED2, Timestamps, info] = load_open_ephys_data('100_ADC6.continuous');

% read the stimulus list text file
StimList = repmat((1:20),1,5);
baseline = 1:4000;

% for each stimulus type
for i = 1:numel(unique(StimList))
    % get indices of all repeats
    reps = find(StimList==StimList(i));
    traces = [];
    % for each repeat
    for j = 1:numel(reps)
        idxON = find(Timestamps==TrialTimeStamps(reps(j),1),1,'first');
        idxOFF = find(Timestamps==TrialTimeStamps(reps(j),2),1,'first');
        nidx = length(idxON:idxOFF);
        mytrace = interp1(1:nidx,PID(idxON:idxOFF),32:32:nidx);
        % subtract baseline
        mytrace = mytrace - mean(mytrace(baseline));
        traces(j,1:length(mytrace)) = mytrace;
    end
    subplot(4,5,i);
    plot(traces');
end

% rescale plots
for i = 1:10
    subplot(4,5,i);
    set(gca,'YLim',[-0.2 0.5]);    
end
for i = 11:15
    subplot(4,5,i);
    set(gca,'YLim',[-0.2 0.6]);    
end
for i = 16:20
    subplot(4,5,i);
    set(gca,'YLim',[-0.2 6]);    
end