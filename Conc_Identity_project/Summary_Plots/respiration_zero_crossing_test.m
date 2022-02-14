%Respiration crossing test

respData_filtered = Respiration_DS_Filt;

% find peak, valley, zci for pressure data

[pks_1,locs_1,w_1,p_1] = findpeaks(respData_filtered, 'MinPeakProminence', 0.3, 'MinPeakDistance', 20);
[pks_2,locs_2,w_2,p_2] = findpeaks(-respData_filtered, 'MinPeakProminence', 0.3, 'MinPeakDistance', 20, 'MinPeakHeight', 0.1);

zci = @(v) find(diff(sign(v))<0 & diff(v) < -0.001);
zero_crossings = zci(respData_filtered);

figure;
plot(respData_filtered);
hold on;
plot(locs_1,respData_filtered(locs_1),'or');
plot(locs_2,respData_filtered(locs_2),'ob');
plot(zero_crossings, respData_filtered(zero_crossings),'og');