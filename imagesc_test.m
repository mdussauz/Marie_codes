% all_trials_response_matrix_of_all_units
% PSTH_matrix_for_all_units

%written by MD 

%% Get Events timestamps from open ephys file
% DON T NEED THAT PART HERE 

% myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with open ephys data
% filename = fullfile(myKsDir,'all_channels.events');
% 
% [Events] = ParseOpenEphysEvents(filename);
% % trial events are in Channel0, odor events are in Channel1
% 
% TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
% TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
% OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];
% 
% TrialTimestamps = TrialTimestamps - offset; %adjust for offset
% OdorTimestamps = OdorTimestamps - offset; % adjust for offset 


%% 
line = 0; % initialize 
step = 0.01;

for whichcluster = 1:length (goodcluster) %*to run multiple clusters
    clear spikecount 
    
  if ~isempty( goodcluster(whichcluster).spikes)  
      %mycluster %checking that it runs through each good clusters
      line = line+1;
      
    for mytrial = 1:NTrials
            firstbin = -10;
            lastbin = 10;
            tbins = firstbin:step:lastbin;
            SpikeTimes = goodcluster(whichcluster).spikes{1, mytrial};
            myspikecount = histcounts(SpikeTimes, tbins);
            spikecount(mytrial,:) = [myspikecount];
    end

    meanalltrialspikecount = ((sum(spikecount, 1)).*(1/step))./NTrials; %sum(A,1) = sum along each column of the matrix A 
    %meanalltrialspikecount %check
    meanalltrial(line,:) = [meanalltrialspikecount];
    

    
    clear meanalltrialspikecount 


  end

    
end 

%size(meanalltrial)
figure(1)
h = imagesc(meanalltrial);
set(h, 'XData', [firstbin, lastbin]);
title('Response matrix for all good units')
xlim([-10, 10])
xlabel('Time')
ylabel('unit #')
colorbar

