% written by MD 

a = readNPY('spike_times.npy');
b = readNPY('spike_clusters.npy');
c = readNPY('amplitudes.npy');

sample_rate = 30000;

which_unit = 5; % to select which unit to check and plot 

spiketimeforunit= double(a(find(b==which_unit))); %get spike timestamps for cluster 5 that I tagged as good
spiketimeforunit= (spiketimeforunit)./ sample_rate; % spike time in s
spikeamplitudefor5= double(c(find(b==5))); %get spike amplitude for cluster 5 that I tagged as good

figure(1)
histogram(spiketimeforunit, 200); %firing rate nb bins = 200
title('FR')

figure(2)
scatter(spiketimeforunit,spikeamplitudefor5, 'filled'); % amplitude across time
title('Amplitude')

%% Read Stimulus File
stimfilename = '190910_17_01.txt';
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = length(unique(StimList));
NTrials= length(StimList);

%% read and process events timestamps from open ephys file
myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with open ephys data
filename = fullfile(myKsDir,'all_channels.events');

% in all_channels there are only timestamps from recorded events 
[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

%Hack for this specific file:
% I found that the offset is - 278.59 => you just look at the first time in
% Timestamps - seems tp agree with value from AdjustClockOffset function 

%offset = 278.59;

% adjust for clock offset between open ephys and kilosort
[offset] = AdjustClockOffset(myKsDir);
offset = offset/sampleRate;
timestamps = timestamps - offset;

TrialTimestamps = TrialTimestamps - offset;
OdorTimestamps = OdorTimestamps - offset; 

[PressureSensor, Timestamps, info] = load_open_ephys_data('100_ADC1.continuous');

Timestamps = Timestamps - offset;

%% raster plot for full session
figure(3)
plot([1 1].*TrialTimestamps(:),[-1 1].*2,'k-');hold on % trial in black
plot([1 1].*OdorTimestamps(:),[-1 1].*2,'r-');hold on %odor in red
plot([1 1].*spiketimeforunit(:),[-0.75 0.75].*2,'b-');hold on

%x =[1 1].*TrialTimestamps(:);
% x1= [1 1 1 1].*TrialTimestamps(:);
% y1= [-1 -1 1 1].*2;
% x2= [1 1].*OdorTimestamps(:);
% y2= [-1 -1 1 1].*2;
% patch(x1(:),y1(:),'EdgeColor','none','FaceColor','grey','FaceAlpha',.3);hold on % trial in black
% patch(x2(:),y2(:),'EdgeColor','none','FaceColor','red','FaceAlpha',.3);hold on %odor in red
% plot([1 1].*spiketimefor5(:),[-0.75 0.75].*2,'b-');hold on


% plotlim = 400; %for checking
% xlim([0 plotlim]) %for checking

xlim([0 max(Timestamps)])
ylim([-5 5])
hold off

%% Assign spikes to trial 

for mytrial = 1:NTrials
    disp(mytrial)
    TrialSpikes = [];
    thisTrialSpikes= spiketimeforunit(find(spiketimeforunit>=TrialTimestamps(mytrial,1) & spiketimeforunit<=TrialTimestamps(mytrial,2)));
    thisOdorSpikes= spiketimeforunit(find(spiketimeforunit>=OdorTimestamps(mytrial,1) & spiketimeforunit<=OdorTimestamps(mytrial,2)));
    %TrialSpikes(mytrial) = thisTrialSpikes;
    
    % PSTH try 
%     figure(4)
%     subplot(4,5,mytrial);
%     plot([1 1].* thisTrialSpikes(:),[-0.75 0.75].*2,'b-')

    for myspike = 1:numel(thisTrialSpikes)
        if ~isempty(thisTrialSpikes)
        plot([1 1].* (thisTrialSpikes(myspike)- TrialTimestamps(mytrial,1)),[0 1] + mytrial ,'k-'); hold on 
        end
    end 

    for myspike2 = 1:numel(thisOdorSpikes)
        if ~isempty(thisOdorSpikes)
        plot([1 1].* (thisOdorSpikes(myspike2) - TrialTimestamps(mytrial,1)),[0 1] + mytrial ,'r-','LineWidth', 1.5); hold on 
        end
    end

end

% plot([OdorTimestamps(1,1) OdorTimestamps(1,1)], [0 NTrials+2], 'r--'); hold on
% plot([OdorTimestamps(1,2) OdorTimestamps(1,1)], [0 NTrials+2], 'r--'); hold on

plot([10 10], [0 NTrials+2], 'r--'); hold on % dotted red line odor on
plot([14 14], [0 NTrials+2], 'r--'); hold on % dotted red line odor off
axis
xlim([0 21])
ylim([0 NTrials+2])
xlabel('Trial Number')
ylabel('Time (s)')
hold off


%% Assign spikes to trial types

% hack to check whether conc (20 stim) or id (16 stim) exp 
if (ismember(20, StimList))==1 %stim 20 should only be present in conc exp
    TrialType=0; %0 for conc exp
    disp('concentration exp')
else 
    TrialType=1; %1 for id exp
    disp('identity exp')
end 

switch TrialType

    case 0 % concentation experiment
        
for mytrial = 1:NTrials
    %disp(mytrial)
    TrialSpikes = [];
    thisTrialSpikes= spiketimeforunit(find(spiketimeforunit>=TrialTimestamps(mytrial,1) & spiketimeforunit<=TrialTimestamps(mytrial,2)));
    thisOdorSpikes= spiketimeforunit(find(spiketimeforunit>=OdorTimestamps(mytrial,1) & spiketimeforunit<=OdorTimestamps(mytrial,2)));
    %TrialSpikes(mytrial) = thisTrialSpikes;
    
    % PSTH try 
%     figure(4)
%     subplot(4,5,mytrial);
%     plot([1 1].* thisTrialSpikes(:),[-0.75 0.75].*2,'b-')


    for myspike = 1:numel(thisTrialSpikes)
        if ~isempty(thisTrialSpikes)
        plot([1 1].* (thisTrialSpikes(myspike)- TrialTimestamps(mytrial,1)),[0 1] + mytrial ,'k-'); hold on 
        end
    end 

    for myspike2 = 1:numel(thisOdorSpikes)
        if ~isempty(thisOdorSpikes)
        plot([1 1].* (thisOdorSpikes(myspike2) - TrialTimestamps(mytrial,1)),[0 1] + mytrial ,'r-','LineWidth', 1.5); hold on 
        end
    end

end

% plot([OdorTimestamps(1,1) OdorTimestamps(1,1)], [0 NTrials+2], 'r--'); hold on
% plot([OdorTimestamps(1,2) OdorTimestamps(1,1)], [0 NTrials+2], 'r--'); hold on

plot([10 10], [0 NTrials+2], 'r--'); hold on % dotted red line odor on
plot([14 14], [0 NTrials+2], 'r--'); hold on % dotted red line odor off
axis
xlim([0 21])
ylim([0 NTrials+2])
xlabel('Trial Number')
ylabel('Time (s)')
hold off

    case 1 % identity experiment 
       