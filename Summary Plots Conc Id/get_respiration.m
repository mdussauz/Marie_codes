dir = '~/Desktop/Analysis/';
experiment = '2021-01-25_09-35-25';


myKsDir = fullfile(dir, experiment);
foo = fullfile(myKsDir, '100_ADC1.continuous');
[Respiration, timestamps, ~] = load_open_ephys_data(foo); % data has channel IDs

% adjust timestamps to account for the start offset in OEPS
OepsSampleRate = 30000; % Open Ephys acquisition rate
%[offset] = AdjustClockOffset(myKsDir);
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

%%
figure;
plot(respData_filtered);
hold on;
plot(locs_1,respData_filtered(locs_1),'or'); %exhalation
plot(locs_2,respData_filtered(locs_2),'ob'); %inhalation
plot(zero_crossings, respData_filtered(zero_crossings),'og'); %approx start of inhalation
