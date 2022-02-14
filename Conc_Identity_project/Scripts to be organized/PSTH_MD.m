%written by MD

mycluster = 1;

% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
offset = offset/sampleRate;
timestamps = timestamps - offset;

%% Get Events 
[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

TrialTimestamps = TrialTimestamps - offset; %adjust for offset
OdorTimestamps = OdorTimestamps - offset; % adjust for offset 

%% PSTH of one unit for one trial 

tstart = TrialTimestamps(1,1); 
tstop = TrialTimestamps(1,2); 

SpikeTimes = goodcluster(mycluster).spikes{1,1};
%SpikeTimes = SpikeTimes.'


step = 0.2; %bins in s
tbins = 0:step:ceil(tstop-tstart); 

spikecount = histcounts(SpikeTimes, tbins); %for cluster 11
spikecount = spikecount .*(1/step);

figure(1)
plot(tbins(1:end-1),spikecount, '-k');  hold on
%histogram(Rate, tbins '-k');
xlim([0 max(tbins)]);

plot([10 10], [0 20], 'r--'); hold on % dotted red line odor on
plot([14 14], [0 20], 'r--'); hold on
hold off

%% PSTH of one unit for all trials


for mytrial = 1:NTrials
	tstart = TrialTimestamps(mytrial,1); 
	tstop = TrialTimestamps(mytrial,2);  
    lastbin = ceil(tstop-tstart); 
    tbins = 0:step:ceil(tstop-tstart); 
    SpikeTimes = goodcluster(mycluster).spikes{1, mytrial};
    myspikecount = histcounts(SpikeTimes, tbins);
    spikecount(mytrial,:) = [myspikecount];
end

meanalltrialspikecount = ((sum(spikecount, 1)).*(1/step))./NTrials; %sum(A,1) = sum along each column of the matrix A 
smoothmeanall = smoothdata(meanalltrialspikecount, 'gaussian', 10);
%test = mean(spikecount); %gets mean of each column 
stdalltrialspikecount = std((sum(spikecount, 1)).*(1/step)) ;
semalltrialspikecount = stdalltrialspikecount ./sqrt(NTrials);

figure(2)
plot(tbins(1:end-1),meanalltrialspikecount, '-k');  hold on
plot(tbins(1:end-1),smoothmeanall, '-c', 'LineWidth', 2);  hold on
%errorbar(semalltrialspikecount)
%histogram(Rate, tbins '-k');
xlim([0 max(tbins)]);
plot([10 10], [0 max(meanalltrialspikecount)+1], 'r--'); hold on % dotted red line odor on
plot([14 14], [0 max(meanalltrialspikecount)+1], 'r--'); hold on
title(['PSTH of cluster #', num2str(goodcluster(mycluster).id), ' for all trials', ' - tbins(s) ',  num2str(step) ])
hold off

%%
figure(3); imagesc(meanalltrialspikecount);

%% PSTH of one unit sorted by concentation and odor identity 
%mycluster = 1;
figure(4)
sgtitle(['PSTH of cluster #', num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step) ])

for i = 1:NOdors %for each odor
    reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
    clear spikecount 

    
    for j = 1:numel(reps)
            
    tstart = TrialTimestamps(reps(j),1); 
	tstop = TrialTimestamps(reps(j),2);  
    lastbin = ceil(tstop-tstart); 
    tbins = 0:step:ceil(tstop-tstart); 
    SpikeTimes = goodcluster(mycluster).spikes{1, reps(j)};
    myspikecount = histcounts(SpikeTimes, tbins);
    spikecount(j,:) = [myspikecount];
    

    end
        
        meanspikecount = ((sum(spikecount, 1)).*(1/step))./numel(reps);
        smoothmean = smoothdata(meanspikecount, 'gaussian', 10);
        subplot(4, 5, StimList(i))
        plot(tbins(1:end-1),meanspikecount, '-k', 'LineWidth', 0.4); hold on
        plot(tbins(1:end-1),smoothmean, '-c', 'LineWidth', 1);  hold on
%         figure(5); subplot(4, 5, StimList(i));  imagesc(meanspikecount);
        
        clear meanspikecount
 end 
    
%totalspike = sum(allspikecount,1)
    for k = 1:20
        
        subplot(4, 5, k)
        
        plot([10 10], [0 40], 'r--'); hold on % dotted red line odor on
        plot([14 14], [0 40], 'r--'); hold on % dotted red line odor off
        
        %axis
        xlim([0 21])
        ylim([0 max(meanalltrialspikecount)+1])
       
    end
    



