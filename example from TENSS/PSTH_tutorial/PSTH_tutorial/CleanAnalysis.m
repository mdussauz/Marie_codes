%clearvars -except TimeMSFromStartTTLs VocTTLON VocOnset  FramesTTLCam2 FrameCountCam2 BehaviorEvent1Frame BehaviorEvent2Frame BehaviorEvent3Frame SpikeTimesInMS


% Hello!
% What do we have? 
% This  workspace has data from an experiment featuring a camera, an
% ultrasound microphone and an ephys system. 

% Camera Stuff:

%  FrameTTLCam2 is the number of frame corresponding to the onset of TTL's
%  FrameCountCam2 is metadata from the camera and corresponds to the FrameNumber of the frames acquired
%  BehaviorEvent1Frame frame numbers for Event 1
%  BehaviorEvent2Frame frame numbers for Event 2
%  BehaviorEvent3Frame frame numbers for Event 3

% Audio Stuff

% VocTTL times in audio Time
% VocOnset vocalization onset times in audio time

% Ephys Stuff

% TimeMSFromStartTTLs   Time for the onset of TTLs in Ephys Time
% SpikeTimesinMS       Time for the spikes in Ephys time


%% Let's get a first glance of our TLL pulses in the different devices

figure; 

subplot(2,2,1)

plot(TimeMSFromStartTTLs, 'ok')

title([ 'Ephys TTL number= ' num2str(numel(TimeMSFromStartTTLs))])



subplot(2,2,2)

plot(FramesTTLCam2, 'or')

title([ 'Camera TTL number= ' num2str(numel(FramesTTLCam2))])


subplot(2,2,3)

plot(VocTTLON, 'og')

title([ 'Audio TTL number= ' num2str(numel(VocTTLON))])


% Ok great, we have the same number of TTLs in all the devices,But they are
% all in different scales
subplot(2,2,4)
hold on
plot(VocTTLON, 'og')
plot(FramesTTLCam2, 'or')
plot(TimeMSFromStartTTLs, 'ok')

%Now we
%can use that to syncronize the Audio and Camera with the Ephys Time

%% First Lets Syncronize the Audio Data with the Ephys Times

figure;

%Lets look at the relationship between TTL's
subplot(2,2,1)

plot(VocTTLON, TimeMSFromStartTTLs, 'og');

xlabel('Audio Time')
ylabel('Ephys Time')
hold on 

plot([VocOnset VocOnset] , ylim, '-k');
%%


VocalizationTimeInEphysTime=interp1(VocTTLON, TimeMSFromStartTTLs, VocOnset, 'linear', 'extrap');


plot(VocOnset , VocalizationTimeInEphysTime, 'ob');

%Great we are N'SYNC, anyone get it or am I too old?

%Lets work with the Camera stuff now.





%% Let's just apply the same Logic for the Camera now


subplot(2,2,2)


plot(FramesTTLCam2, TimeMSFromStartTTLs, 'or');
xlabel('Frame Number')
ylabel('Ephys Time')

hold on


plot([FrameCountCam2(BehaviorEvent1Frame) FrameCountCam2(BehaviorEvent1Frame)] , ylim, '-k');
plot([FrameCountCam2(BehaviorEvent2Frame) FrameCountCam2(BehaviorEvent2Frame)] , ylim, '-c');
plot([FrameCountCam2(BehaviorEvent3Frame) FrameCountCam2(BehaviorEvent3Frame)] , ylim, '-m');



FrameTimeinEphys=interp1(FramesTTLCam2, TimeMSFromStartTTLs, FrameCountCam2, 'linear', 'extrap');



Event1TimeSync=FrameTimeinEphys(BehaviorEvent1Frame);
Event2TimeSync=FrameTimeinEphys(BehaviorEvent2Frame);
Event3TimeSync=FrameTimeinEphys(BehaviorEvent3Frame);

plot(FrameCountCam2(BehaviorEvent1Frame), Event1TimeSync, 'ok');
plot(FrameCountCam2(BehaviorEvent2Frame), Event2TimeSync, 'oc');
plot(FrameCountCam2(BehaviorEvent3Frame), Event3TimeSync, 'om');


%% Great! All behavioral variables are syncronized to the Ephys! Now lets see what the cell is doing






figure;

% Lets fist just bin the Spikes for the cell in 1000ms bins
step=1000;
tbins = FrameTimeinEphys(1):step :FrameTimeinEphys(end); % time bin centers for spike train binnning
Rate = hist(SpikeTimesInMS,tbins); %Histogram for rate in time

subplot(2,2,1:2)
plot(tbins,Rate, '-k');
xlim([0 max(tbins)]);
ylabel('Rate (Hz)')
xlabel('Time (s)')
hold on 

% LEts plot on top our behavioral events
plot([Event1TimeSync Event1TimeSync], ylim, '-k')
plot([Event2TimeSync Event2TimeSync], ylim, '-c')
plot([Event3TimeSync Event3TimeSync], ylim, '-m')

%% And maybe now the vocalizations
plot([VocalizationTimeInEphysTime VocalizationTimeInEphysTime], ylim, '-g')


%%  

subplot(2,2,3:4)

imagesc(Rate)
hold on 

plot([Event1TimeSync/1000 Event1TimeSync/1000], ylim, '-w')
plot([Event2TimeSync/1000 Event2TimeSync/1000], ylim, '-c')
plot([Event3TimeSync/1000 Event3TimeSync/1000], ylim, '-m')


% The cells seems to be activated after Event 1

%% Let's do a PSTH

      bw=50;    %in ms
      PreTime=3;  %in sec
      PostTime=3;  %in sec
      
      % compute plotting times
lastBin = ceil((PostTime*1000)); % last bin edge in ms
edge=-PreTime*1000:bw:lastBin; %
xmin = edge(1); % for plotting
xmax = edge(end); % for plotting
ntrials=length(Event3TimeSync);

      psth = zeros(numel(edge),length(Event3TimeSync));
       for trial = 1:length(Event3TimeSync)
           psthSpikes =SpikeTimesInMS(((Event3TimeSync(trial)/1000-PreTime)) < SpikeTimesInMS/1000  &  SpikeTimesInMS/1000 < (Event3TimeSync(trial)/1000 +PostTime)); % extract spikes occuring within pre/post times of start time
           psth(:,trial) = histc((psthSpikes-Event3TimeSync(trial)),edge) / (bw/1000); %/ ntrials; % make spikes relative to start and extract spike count per bin, divide by bw in seconds to get firing rate 
           clear psthSpikes
       end
       psthTrialAvg = mean(psth,2); % take mean across trials
       psthylim = [0 max(psthTrialAvg)]; % for plotting
       
       
       figure;
        subplot(2,1,1);
       bar(edge,psthTrialAvg,'histc');
       
       xlim([xmin xmax])
       ylim(psthylim)  
       title(' Event3 PSTH for cell');
       ylabel('Spike Rate');

       
