function  BasicCheckStim(myKsDir)
% written by MD

% Do basic checks on recording session and associated stim file: 
% - whether session is an identity or concentration experiment
% - check that the number of trials during the session match the settings
% - check that the duration of trial ON and odor ON match the settings
% Figure 1: plot for each odor when they were on during session
% Figure 2: plot trial ON/OFF odor ON/OFF and odor number

% as we are not looking at kilosort data, there is no need to adjust the
% clock offset between open ephys and kilosort data) 


%clc; clear all; close all;

%% add the relevant folders and repositories to path   
addpath(genpath('/mnt/data')) % path to data folder 
addpath(genpath('/opt/open-ephys-analysis-tools'))% path to open ephys scripts
%addpath(genpath('/home/marie/Documents/Open Ephys code/analysis-tools-master'))% path to open ephys scripts

%% Filepaths
myKsDir = '/mnt/data/N5/2019-09-10_17-01-25'; % directory with open ephys data

%myKsDir = '/mnt/data/N6/2019-09-09_14-02-33'

%% read and process events timestamps from open ephys file
filename = fullfile(myKsDir,'all_channels.events');
% in all_channels there are only timestamps from recorded events 
[Events] = ParseOpenEphysEvents(filename);
% trial events are in Channel0, odor events are in Channel1

TrialTimestamps = [Events.Channel0.On Events.Channel0.Off];
TrialTimestamps(1,:) = []; % delete first entry - just an empty trigger
OdorTimestamps = [Events.Channel1.On Events.Channel1.Off];

%% read data from open ephys file - get timestamps
%datafilename = input('Enter data folder name: ' , 's');
[PressureSensor, Timestamps, info] = load_open_ephys_data('100_ADC1.continuous');

%% read the stimulus txt file 
stimfilename = input('Enter Stimulus file name: ' , 's'); % eg 190910_17_01.txt !no quotes
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);

%Info from txt file
NOdors = numel(unique(StimList)); %nber of odors

%% hack to check whether conc (20 stim) or id (16 stim) exp 
if (ismember(20, StimList))==1 %stim 20 should only be present in conc exp
    TrialType=0; %0 for conc exp
    disp('concentration exp')
else 
    TrialType=1; %1 for id exp
    disp('identity exp')
end 

%% get number of trials and check with settings
NTrials = (numel(OdorTimestamps))./2; %actual value 
setNTrials= numel(StimList); %from settings
if NTrials==setNTrials
    disp('Nber of trials is correct')
else 
    disp('Nber of trials is incorrect')
end


%% get Trial mean Time and Odor mean time  and check with settings
meanTrialTime = mean(TrialTimestamps(:,2)-TrialTimestamps(:,1)); % actual value %col1 is on and col2 is off
meanOdorTime = mean(OdorTimestamps(:,2)-OdorTimestamps(:,1)); % actual value
setTrialTime = (StimTime(1)+StimTime(2)+StimTime(3))./1000; % sec from settings
setOdorTime = StimTime(2)./1000; % sec from settings


if fix(meanTrialTime)==setTrialTime && fix(meanOdorTime)==setOdorTime
    disp('Trial and odor times are correct')
else 
    disp('Trial and/or odor times are incorrect')
end

%% Assign odor number to each trial 
for i = 1:NOdors %for each odor
    reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
    %SortedTrials = [];
    %initializing
    idxON = zeros(1,numel(reps)); 
    idxOFF = zeros(1,numel(reps)); 
    OdidxON = zeros(1,numel(reps)); 
    OdidxOFF = zeros(1,numel(reps)); 
    

    for j = 1:numel(reps) % for each repeat %more flexible to errors this way than using Nrepeats
        %get index of timestamps for each trial ON/OFF
        idxON(j) = find(Timestamps==TrialTimestamps(reps(j),1),1,'first');
        idxOFF(j) = find(Timestamps==TrialTimestamps(reps(j),2),1,'first');
        
        %get index of timestamps for each odor ON/OFF
        OdidxON(j) = find(Timestamps==OdorTimestamps(reps(j),1),1,'first');
        OdidxOFF(j) = find(Timestamps==OdorTimestamps(reps(j),2),1,'first');
    end
    
    %store for each odor number the index of the Timestamps at which the
    %trial was on
    SortedTrial.(['Odor',num2str(StimList(i))]).On =idxON;
    SortedTrial.(['Odor',num2str(StimList(i))]).Off =idxOFF;
    
        
    %store for each odor number the index of the Timestamps at which they
    %happened
    SortedOdor.(['Odor',num2str(StimList(i))]).On = OdidxON;
    SortedOdor.(['Odor',num2str(StimList(i))]).Off =OdidxOFF;
    
    %store for each odor number the trial number of each of their repeat
    repeats.(['Odor',num2str(StimList(i))])=reps;
    
    %making different plots for each odor - playing with graphs 
     figure(1)
     subplot(NOdors,1,i);
     
     %using the repeat number
     plot([1 1].*TrialTimestamps(repeats.(['Odor',num2str(StimList(i))]),1),[-1 1].*2,'b-');hold on %trial in blue
     plot([1 1].*TrialTimestamps(repeats.(['Odor',num2str(StimList(i))]),2),[-1 1].*2,'b-');hold on %odor in blue
     plot([1 1].*OdorTimestamps(repeats.(['Odor',num2str(StimList(i))]),1),[-1 1].*2,'r-');hold on %odor in red
     plot([1 1].*OdorTimestamps(repeats.(['Odor',num2str(StimList(i))]),2),[-1 1].*2,'r-');hold off %odor in red
     
     %following lines are a different way to do the same thing but
     %using Sorted structure instead
%      plot([1 1].*Timestamps(SortedTrial.(['Odor',num2str(StimList(i))]).On),[-1 1].*2,'b-');hold on % on - blue for trial with spe odor  
%      plot([1 1].*Timestamps(SortedTrial.(['Odor',num2str(StimList(i))]).Off),[-1 1].*2,'b-');hold on % off - blue for trial with spe odor
%      plot([1 1].*Timestamps(SortedOdor.(['Odor',num2str(StimList(i))]).On),[-1 1].*2,'r-');hold on % on - blue for trial with spe odor  
%      plot([1 1].*Timestamps(SortedOdor.(['Odor',num2str(StimList(i))]).Off),[-1 1].*2,'r-');hold off % off - blue for trial with spe odor
 
      
     
     xlim([0 max(Timestamps)])
     ylim([-4 4])
     ylabel(['Od',num2str(StimList(i))])
    
  
end

%% Plot trial on/off and odor on/off

figure(2)
plot([1 1].*TrialTimestamps(:),[-1 1].*2,'k-');hold on % trial in black
plot([1 1].*OdorTimestamps(:),[-1 1].*2,'r-');hold on %odor in red
xlim([0 max(Timestamps)])
ylim([-4 4])
%annotating odor number on top - playing around with graphs and checking
%code
for l = 1:NOdors
    for m = 1:numel(reps)
text(Timestamps(SortedTrial.(['Odor',num2str(StimList(l))]).On(m)),2.5,['Od',num2str(StimList(l))])
    end
end
hold off


%% Reorganize conc odor nber 
%in progress

end 