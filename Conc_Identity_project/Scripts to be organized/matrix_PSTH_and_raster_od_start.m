%matrix_PSTH_and_raster_od_start

%% Read stim file 

stimfilename = '190918_16_35.txt';


[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList);

%% PSTH of one unit sorted by concentation and odor identity 
mycluster = 4;
figure(1)
sgtitle(['PSTH of cluster #', num2str(goodcluster(mycluster).id), ' sorted by trial type', ' - tbins(s) ',  num2str(step) ])

for i = 1:NOdors %for each odor
    reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
    clear spikecount 

    
    for j = 1:numel(reps)
               
    firstbin = -10;
    lastbin = 10;
    tbins = firstbin:step:lastbin;
    SpikeTimes = goodcluster(mycluster).spikes{1, reps(j)};
    myspikecount = histcounts(SpikeTimes, tbins);
    spikecount(j,:) = [myspikecount];

    end
        
        meanspikecount = ((sum(spikecount, 1)).*(1/step))./numel(reps);
        
        taxis = -500:500;  % make a time axis of 1000 ms
        t_wid = 100;  % width of kernel (in ms)
        gauss_kernel = normpdf(taxis, 0, t_wid);
        gauss_kernel = gauss_kernel ./ sum(gauss_kernel);
        zs = conv(meanspikecount,gauss_kernel,'same');
        %FR = zs*1000; %converting firing rate to Hz.
        FR = zs;
        % FR contains smoothed firing rates;
        
        subplot(4, 5, StimList(i))
        %plot(tbins(1:end-1),meanspikecount, '-c', 'LineWidth', 0.4); hold on
        plot(tbins(1:end-1),FR, '-k', 'LineWidth', 0.6);  hold on
        
        clear meanspikecount
 end 
    
%totalspike = sum(allspikecount,1)
    for k = 1:20
        
        subplot(4, 5, k)
        
        plot([0 0], [0 125], 'r--'); hold on % dotted red line odor on
        plot([4 4], [0 125], 'r--'); hold on % dotted red line odor off
        
        %axis
        xlim([-10 10])
        ylim([0 125])
       
    end
    

subplot(4, 5, 1)
title('Odor 1', 'FontSize',8 ) 
ylabel('Conc 10^-4', 'FontSize',8)
subplot(4, 5, 2)
title('Odor 2', 'FontSize',8 ) 
subplot(4, 5, 3)
title('Odor 3', 'FontSize',8 ) 
subplot(4, 5, 4)
title('Odor 4', 'FontSize',8 ) 
subplot(4, 5, 5)
title('Odor 5', 'FontSize',8 ) 
subplot(4, 5, 6)
ylabel('Conc 10^-3', 'FontSize',8)
subplot(4, 5, 11)
ylabel('Conc 10^-2', 'FontSize',8)
subplot(4, 5, 16)
xlabel('Time (s)', 'FontSize',8)
ylabel({'Conc 10^-1','Repeat #'}, 'FontSize',8)

hold off





% %% Plot raster matrix sorted by concentation and odor identity 
% % matrix of raster plots sorted by odor type on x axis and concentration on
% 
% 
% % y axis
% 
% % color settings to check repeats
% %color = {'r-','m-','b-','g-', 'y'}
% 
% %for mycluster = 1:length (goodcluster) *to run multiple clusters
%  mycluster = 4; %which cluster to plot
%  
%     if ~isempty( goodcluster(mycluster).spikes)
%     %figure(mycluster)     
%     figure(2)
%     sgtitle(['cluster #', num2str(goodcluster(mycluster).id)])   
%     
% 
%         
%     for i = 1:NOdors %for each odor
%     reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
%     
%     for j = 1:numel(reps)
%         
%        
%             
%              spikearray = gpuArray(goodcluster(mycluster).spikes{1,reps(j)}); %if gpu
%            %for spike = 1:numel(goodcluster(mycluster).spikes{1,reps(j)}) %use this line instead if no gpu:
%    
%             
%                 if ~isempty( goodcluster(mycluster).spikes{1, reps(j)})
%                     subplot(4, 5, StimList(i))
%                     %plot([1 1].* (goodcluster(mycluster).spikes{1,reps(j)}),[0 1] + j-1 ,color{j}); hold on %to check repeats
%                     plot([1 1].* (goodcluster(mycluster).spikes{1,reps(j)}),[0 1] + j-1 ,'k-'); hold on %if no gpu
%                    plot([1 1].* (spikearray),[0 1] + j-1,'k-','LineWidth', 0.5); hold on
%                    
%                    %plot(spikearray, (j-1), 'o')
%                 end 
%             
%             
% 
%         end 
%     end 
%      
%  disp(['Cluster ', num2str(mycluster), ' is done'])      
%     else 
%         disp(['entry ', num2str(mycluster), ' is empty'])
%     end
%     
%     for l = 1:20
%         subplot(4, 5, l)
%         
%         plot([0 0], [0 6], 'r--'); hold on % dotted red line odor on
%         plot([4 4], [0 6], 'r--'); hold on % dotted red line odor off
% %         plot([10 10], [0 6], 'k--'); hold on % dotted black line odor on
% %         plot([14 14], [0 6], 'k--'); hold on % dotted black line odor off
%         
%         %axis
%         xlim([-10 10])
%         ylim([0 5])
%        
%     end
%     
% subplot(4, 5, 1)
% title('Odor 1', 'FontSize',8 ) 
% ylabel('Conc 10^-4', 'FontSize',8)
% subplot(4, 5, 2)
% title('Odor 2', 'FontSize',8 ) 
% subplot(4, 5, 3)
% title('Odor 3', 'FontSize',8 ) 
% subplot(4, 5, 4)
% title('Odor 4', 'FontSize',8 ) 
% subplot(4, 5, 5)
% title('Odor 5', 'FontSize',8 ) 
% subplot(4, 5, 6)
% ylabel('Conc 10^-3', 'FontSize',8)
% subplot(4, 5, 11)
% ylabel('Conc 10^-2', 'FontSize',8)
% subplot(4, 5, 16)
% xlabel('Time (s)', 'FontSize',8)
% ylabel({'Conc 10^-1','Repeat #'}, 'FontSize',8)
% 
% 
%  hold off 
%  
% 
% %end *to run multiple clusters

