%raster_2dmatrix_sorted

%written by MD

% matrix of raster plots sorted by odor type on x axis and concentration on
% y axis

% color settings to check repeats
%color = {'r-','m-','b-','g-', 'y'}

%% Read stim file 
%session 1 - conc
stimfilename = '190910_17_01.txt';

[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList);

%% Plot matrix

%for mycluster = 1:length (goodcluster) *to run multiple clusters
 mycluster = 1; %which cluster to plot
 
    if ~isempty( goodcluster(mycluster).spikes)
    figure(mycluster)     
    sgtitle(['cluster #', num2str(goodcluster(mycluster).id)])   
    

        
    for i = 1:NOdors %for each odor
    reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
    
    for j = 1:numel(reps)
        
       
            
            for spikearray = gpuArray(goodcluster(mycluster).spikes{1,reps(j)}) %if gpu
           %for spike = 1:numel(goodcluster(mycluster).spikes{1,reps(j)}) %use this line instead if no gpu:
   
            
                if ~isempty( goodcluster(mycluster).spikes{1, reps(j)})
                    subplot(4, 5, StimList(i))
                    %plot([1 1].* (goodcluster(mycluster).spikes{1,reps(j)}),[0 1] + j-1 ,color{j}); hold on %to check repeats
                    %plot([1 1].* (goodcluster(mycluster).spikes{1,reps(j)}),[0 1] + j-1 ,'k-'); hold on %if no gpu
                    plot([1 1].* (spikearray),[0 1] + j-1 ,'k-'); hold on 
                 
                end 
            
            end

        end 
    end 
     
 disp(['Cluster ', num2str(mycluster), ' is done'])      
    else 
        disp(['entry ', num2str(mycluster), ' is empty'])
    end
    
    for k = 1:20
        subplot(4, 5, k)
        
        plot([10 10], [0 6], 'r--'); hold on % dotted red line odor on
        plot([14 14], [0 6], 'r--'); hold on % dotted red line odor off
%         plot([10 10], [0 6], 'k--'); hold on % dotted black line odor on
%         plot([14 14], [0 6], 'k--'); hold on % dotted black line odor off
        
        %axis
        xlim([0 21])
        ylim([0 5])
       
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
 

%end *to run multiple clusters