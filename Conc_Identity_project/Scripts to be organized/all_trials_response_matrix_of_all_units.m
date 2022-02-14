% all_trials_response_matrix_of_all_units
% PSTH_matrix_for_all_units

%written by MD 

%% Get Events timestamps from open ephys file
myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with open ephys data
filename = fullfile(myKsDir,'all_channels.events');

[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

TrialTimestamps = TrialTimestamps - offset; %adjust for offset
OdorTimestamps = OdorTimestamps - offset; % adjust for offset 

%% 
line = 0; % initialize 
step = 0.2;

for mycluster = 1:length (goodcluster) %*to run multiple clusters
    clear spikecount 
    
  if ~isempty( goodcluster(mycluster).spikes)  
      %mycluster %checking that it runs through each good clusters
      line = line+1;
      
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
    %meanalltrialspikecount %check
    meanalltrial(line,:) = [meanalltrialspikecount];
    
    %stdalltrialspikecount = std((sum(spikecount, 1)).*(1/step)) ;
     size(spikecount(:))
    z_score = (spikecount(:) - mean(spikecount(:)))./  std((spikecount(:)), 1);
    z_scoreall(line,:) = [z_score];
    
    clear meanalltrialspikecount 
    clear z_score

  end

    
end 

%size(meanalltrial)
figure(1)
h = imagesc(meanalltrial);
set(h, 'XData', [0, lastbin]);
title('Response matrix for all good units')
xlim([0, 20.2])
xlabel('Time')
ylabel('unit #')
colorbar

%size(z_scoreall)
figure(2)  % z_score needs to be fixed 
smoothzscore = smoothdata(z_scoreall, 'gaussian', 8);
g = imagesc(z_scoreall);
%g = imagesc(smoothzscore);
set(g, 'XData', [0, lastbin]);
title('Z score matrix for all good units')
xlim([0, 20.2])
xlabel('Time')
ylabel('unit #')
caxis([-1 1])
colorbar

clear meanalltrial
clear z_scoreall 