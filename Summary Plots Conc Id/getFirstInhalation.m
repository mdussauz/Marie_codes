function [closest_neighbor_zci] = getFirstInhalation(myKsDir, OdorTimestamps)

% dir = '~/Desktop/Analysis/';
% experiment = '2021-01-25_09-35-25';
% myKsDir = fullfile(dir, experiment);

foo = fullfile(myKsDir, '100_ADC1.continuous');
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs

% adjust timestamps to account for the start offset in OEPS
offset = timestamps(1);
timestamps = timestamps - offset;
%%
% Get respiration
newTimeBin = 0.001; % in seconds
newTimeBase = (0:newTimeBin:max(timestamps))';
RespData = interp1q(timestamps,Respiration,newTimeBase);

% rescale the data
RespData = RespData - median(RespData);

% invert: inhalations should be negative
RespData = -RespData;

% smooth pressure data by moving mean filter
respData_filtered = smoothdata(RespData);

% find peak, valley, zci for pressure data

[pks_1,locs_1,w_1,p_1] = findpeaks(respData_filtered, 'MinPeakProminence', 0.3, 'MinPeakDistance', 20);
[pks_2,locs_2,w_2,p_2] = findpeaks(-respData_filtered, 'MinPeakProminence', 0.3, 'MinPeakDistance', 20, 'MinPeakHeight', 0.1);

zci = @(v) find(diff(sign(v))<0 & diff(v) < -0.001);
zero_crossings = zci(respData_filtered);

%% Find nearest zero crossing time right after odor start
odorStart = OdorTimestamps(:,1);
odorStart = odorStart - offset; % bc newTimeBase accounts for offset
zci_times = newTimeBase(zero_crossings);

closest_neighbor_zci = zeros(1, length(odorStart));
for i = 1:length(odorStart)
    dist = zeros(1, length(zci_times));
    for j = 1:length(zci_times)
        dist(j) = zci_times(j) - odorStart(i);
    end
    dist(dist <0) = inf;
    [~,closestIndex] = min(dist);
    closest_neighbor_zci(i) = zci_times(closestIndex);

end

end