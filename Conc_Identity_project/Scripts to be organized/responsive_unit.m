%responsive_unit
%written by MD 

%% Read Stimulus File
stimfilename = '190910_17_01.txt';
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);
NOdors = length(unique(StimList));
NTrials= length(StimList);

%% Get Events timestamps from open ephys file
myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with open ephys data
filename = fullfile(myKsDir,'all_channels.events');

[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

% adjust for clock offset between open ephys and kilosort
sampleRate = 30000; 
[offset] = AdjustClockOffset(myKsDir);
offset = offset/sampleRate;
%timestamps = timestamps - offset;

TrialTimestamps = TrialTimestamps - offset; %adjust for offset
OdorTimestamps = OdorTimestamps - offset; % adjust for offset 

%%  1st take 10 s pre stim and cut in 4s chunks 


mycluster =11;  % choosing cluster
%mytrial = 1; % choosing 1 trial 
iteration = 100;

for mytrial = 1:NTrials
    SpikeTimes = goodcluster(mycluster).spikes{1,mytrial};

    tstart = TrialTimestamps(mytrial,1); 
    tstop = TrialTimestamps(mytrial,2);  
    odorstart = OdorTimestamps(mytrial,1);
    odorstop = OdorTimestamps(mytrial,2);


    prestim = odorstart - tstart; %full duration of prestim period
    timewindow = 4; % 4s

    a = tstart;
    b = odorstart - timewindow; %arbitrary - as this will be sampled via +4 after
    randtstart = (b-a).*rand(iteration,1) + a; % generates 100 integers between tstart and odorstart - 4s
    randtstop = randtstart + 4;

    %align to tstart
    randtstart = randtstart - tstart;
    randtstop = randtstop - tstart; 
    
    odorstart = odorstart - tstart;
    odorstop = odorstop - tstart;

    %hist(randtstart)
    %spikecount = zeros(100, 1000);
    %allFR = zeros(100,1);

    %% constructing distribution of FR in period of 4s for prestim  
    for i = 1:iteration
        clear myspikecount
        %i
        prestimspike = SpikeTimes(find(SpikeTimes>=randtstart(i) & SpikeTimes<=randtstop(i)));

        tbins = randtstart(i):1:randtstop(i);
        myspikecount = histcounts(prestimspike, tbins);
        FR = mean(myspikecount);
        allFR_prestim (mytrial,i) = FR;
        
        %myspikecount = myspikecount .*(1/step);
        %spikecount(mytrial,4) =  [myspikecount];
        %allFR(i) = FR;
    end 

    
    %% getting FR for stim
    
    stimspike = SpikeTimes(find(SpikeTimes>=odorstart & SpikeTimes<=odorstop));
    tbins_odor = odorstart:1:odorstop;
    myspikecount_stim = histcounts(stimspike, tbins_odor);
    FR_stim = mean(myspikecount_stim);
    allFR_stim (mytrial) = FR_stim;
   
end 

figure(1)
x = histogram(allFR_prestim);



figure(2)
y = histogram(allFR_stim);

%% Rank sum test 

p = ranksum(x.Values ,y.Values);
p

if p < 0.05
    disp (['Unit ', num2str(goodcluster(mycluster).id), ' is responsive'])
else 
    disp (['Unit ', num2str(goodcluster(mycluster).id), ' is not responsive'])
end 


%% 