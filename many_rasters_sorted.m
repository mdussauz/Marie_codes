%many_rasters_sorted

% written by MD
%produces 2 matrices of rasters for 1 or many units 

% output = 2 figures:

%Figure 1 = raster across concentrations
% each sublplot = 1 concentration
% x axis = time
% y axis = all repeats all odors - every 5 lines = 1 repeat

%Figure 2 = raster across odor identities
% each sublplot = 1 odor identity 
% x axis = time
% y axis = all repeats all odors - every 5 lines = 1 repeat




%% raster plot sorted by concentration
%for mycluster = 1:length (goodcluster) *to run multiple clusters
 mycluster = 1; %just for checking purpose     
 repconc1 = 0; 
 repconc2 = 0; 
 repconc3 = 0; 
 repconc4 = 0; 
    if ~isempty( goodcluster(mycluster).spikes)
    figure(mycluster);     
    sgtitle(['cluster #', num2str(goodcluster(mycluster).id)])    
        
        for mytrial = 1:NTrials
            
            if ismember(StimList(mytrial), Trial.Conc(1,:))
                conc = 1;
                repconc = repconc1;
                repconc1 = repconc1 + 1;
            elseif ismember(StimList(mytrial), Trial.Conc(2,:))
                conc = 2;
                repconc = repconc2;
                repconc2 =  repconc2 + 1;
            elseif ismember(StimList(mytrial), Trial.Conc(3,:))
                conc = 3;
                repconc = repconc3; 
                repconc3 = repconc3 + 1;
            elseif ismember(StimList(mytrial), Trial.Conc(4,:))
                conc = 4;
                repconc = repconc4;
                repconc4 = repconc4 + 1;
            end
            
            for spike = gpuArray(goodcluster(mycluster).spikes{1, mytrial}))
            %for spike = 1:numel(goodcluster(mycluster).spikes{1, mytrial}) %use this if no gpu 
           
            
                if ~isempty( goodcluster(mycluster).spikes{1, mytrial})
                    subplot(2, 2, conc)
                    plot([1 1].* (goodcluster(mycluster).spikes{1, mytrial}(spike)),[0 1] + repconc ,'k-'); hold on 
                 
                end 
            
            end

        end 

  
 disp(['Cluster ', num2str(mycluster), ' is done'])      
    else 
        disp(['entry ', num2str(mycluster), ' is empty'])
    end
    
    for i = 1:4
        subplot(2, 2, i)
        
        plot([10 10], [0 25], 'r--'); hold on % dotted red line odor on
        plot([14 14], [0 25], 'r--'); hold on % dotted red line odor off
        
        %axis
        xlim([0 21])
        ylim([0 25])
        xlabel('Time (s)')
        ylabel('Repeat Number')
    end
    
subplot(2, 2, 1)
title('conc 10^-4') 
subplot(2, 2, 2)
title('conc 10^-3') 
subplot(2, 2, 3)
title('conc 10^-2') 
subplot(2, 2, 4)
title('conc 10^-1') 
 hold off 
 
 %% raster plot sorted by identity
%for mycluster = 1:length (goodcluster)
 mycluster = 1; %just for checking purpose     
 repid1 = 0; 
 repid2 = 0; 
 repid3 = 0; 
 repid4 = 0; 
 repid5 = 0;
 
    if ~isempty( goodcluster(mycluster).spikes)
    figure(mycluster + 1);     
    sgtitle(['cluster #', num2str(goodcluster(mycluster).id)])    
        
        for mytrial2 = 1:NTrials
            
            if ismember(StimList(mytrial2), Trial.Id(1,:))
                idy = 1;
                repid = repid1;
                repid1 = repid1 + 1;
            elseif ismember(StimList(mytrial2), Trial.Id(2,:))
                idy = 2;
                repid = repid2;
                repid2 =  repid2 + 1;
            elseif ismember(StimList(mytrial2), Trial.Id(3,:))
                idy = 3;
                repid = repid3; 
                repid3 = repid3 + 1;
            elseif ismember(StimList(mytrial2), Trial.Id(4,:))
                idy = 4;
                repid = repid4;
                repid4 = repid4 + 1;
            elseif ismember(StimList(mytrial2), Trial.Id(5,:))
                idy = 5;
                repid = repid5;
                repid5 = repid5 + 1;
            end
            
            for spike2 = 1:numel(goodcluster(mycluster).spikes{1, mytrial2})
           
            
                if ~isempty( goodcluster(mycluster).spikes{1, mytrial2})
                    subplot(3, 2, idy)
                    plot([1 1].* (goodcluster(mycluster).spikes{1, mytrial2}(spike2)),[0 1] + repid ,'k-'); hold on 
                 
                end 
            
            end

        end 

  
 disp(['Cluster ', num2str(mycluster), ' is done'])      
    else 
        disp(['entry ', num2str(mycluster), ' is empty'])
    end
    
    for j = 1:5
        subplot(3, 2, j)
        title(['id', num2str(j)])   
        plot([10 10], [0 20], 'r--'); hold on % dotted red line odor on
        plot([14 14], [0 20], 'r--'); hold on % dotted red line odor off
        
        %axis
        xlim([0 21])
        ylim([0 20])
        xlabel('Time (s)')
        ylabel('Repeat Number')
    end
 hold off 
%end *to run multiple clusters