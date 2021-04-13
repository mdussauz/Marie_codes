% PSTH of one unit sorted by concentation and odor identity  - NORMALIZED

%remove z_score and ormalize by baseline instead 

mycluster = 1;



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
        meanspikecount = mean(spikecount, 1);
        size(meanspikecount)
        stdspikecount = std(spikecount(:),1);
        size(stdspikecount)
        z_score = ((sum(spikecount, 1)) -  meanspikecount )./  stdspikecount;
        size(z_score)
        %smoothmean = smoothdata(meanspikecount, 'gaussian', 10);
        subplot(4, 5, StimList(i))
        %plot(tbins(1:end-1),meanspikecount, '-k', 'LineWidth', 0.4); hold on
        %plot(tbins(1:end-1),smoothmean, '-c', 'LineWidth', 1);  hold on
        plot(tbins(1:end-1),z_score, '-k', 'LineWidth', 1);  hold on
%       figure(5); subplot(4, 5, StimList(i));  imagesc(meanspikecount);
        
        clear meanspikecount
        clear stdspikecount
        clear z_score
 end 
    
%totalspike = sum(allspikecount,1)
    for k = 1:20
        
        subplot(4, 5, k)
        
        plot([10 10], [-10 10], 'r--'); hold on % dotted red line odor on
        plot([14 14], [-10 10], 'r--'); hold on % dotted red line odor off
        
        %axis
        xlim([0 21])
        ylim([-10 10])
       
    end
    
