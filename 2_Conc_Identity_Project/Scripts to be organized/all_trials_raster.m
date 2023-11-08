% all_trials_raster

% written by MD

% plot spikes for each units during all tials 
% y axis trial 1 to last trial 
% x axis time 

%% Read stim file 
%session 1 - conc
stimfilename = '190910_17_01.txt';

[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = numel(unique(StimList));
NTrials= numel(StimList);


%% Plot spikes during trial 

%i = 0; %initialize 

for mycluster = 1:length (goodcluster) %to run multiple clusters
%mycluster = 9;   
clear spikearray
   
    if ~isempty( goodcluster(mycluster).spikes)
        %i = 1+1;   
        figure(mycluster); %mycluster
        plot([10 10], [0 NTrials+2], 'r--'); hold on % dotted red line odor on
        plot([14 14], [0 NTrials+2], 'r--'); hold on % dotted red line odor off
        title(['#', num2str(goodcluster(mycluster).id)])
        %axis
        xlim([0 21])
        ylim([0 NTrials+2])
        xlabel('Time (s)')
        ylabel('Trial Number')
        
        for mytrial = 1:NTrials
            
            spikearray = gpuArray(goodcluster(mycluster).spikes{1, mytrial}); % if gpu
            
            for spike = 1:numel( spikearray) % if gpu
            %for spike = 1:numel(goodcluster(mycluster).spikes{1, mytrial})% use this if no gpu
            
           
            
                %if ~isempty( goodcluster(mycluster).spikes{1, mytrial}) % use this if no gpu
                if ~isempty( spikearray) % use if  gpu
                    
                    %plot([1 1].* (goodcluster(mycluster).spikes{1, mytrial}(spike)),[0 1] + mytrial ,'k-'); hold on % use this if no gpu
                    plot([1 1].* (spikearray(spike)),[0 1] + mytrial ,'k-'); hold on  % if gpu
                 
                end 
            
            end

        end 

 hold off   
 disp(['Cluster', num2str(mycluster), ' is done'])      
    else 
        disp(['entry', num2str(mycluster), 'is empty'])
    end
    

end %to run multiple clusters




