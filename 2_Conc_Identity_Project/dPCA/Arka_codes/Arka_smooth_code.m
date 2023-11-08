taxis = -500:500;  % make a time axis of 1000 ms
t_wid = 100;  % width of kernel (in ms)
gauss_kernel = normpdf(taxis, 0, t_wid);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

clusterNum = 1:29;

for c = 1:length(clusterNum) % loop through each unit
    clusterIdx = clusterNum(c);
    for t = 1:size(times,1) % loop through trials
        tempPSTH = squeeze(PSTH(clusterIdx,t,:));

        %PSTH is where I store spike times (1 for spikes and 0 elsewhere)
        %PSTH is Neurons X trials X time (in ms)

        zs = conv(tempPSTH,gauss_kernel,'same');
        FR(clusterIdx,t,:) = zs*1000; %converting firing rate to Hz.

        % FR contains smoothed firing rates;

    end
end
----------