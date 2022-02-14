function [TrialType, odor_repeats, Trial] = SortTrialType(stimfilename)

%written by MD

%Output:
% - TrialType: 
%Specify what was the experiment 
% 0 for conc exp and 1 for  16 odors exp
% - odor_repeats
% gives the indices of all repeats for each odor 
% - Trial 
%Classify stim number into different odor id and conc based on experiment type

%% get info from stimulus file 
stimfilename = '190910_17_01.txt'; 
[StimTime, StimList, Nrepeats] = ReadStimFile(stimfilename);

NOdors = numel(unique(StimList)); %nber of odors

%% hack to check whether conc (20 stim) or id (16 stim) exp 
if (ismember(20, StimList))==1 %stim 20 should only be present in conc exp
    TrialType=0; %0 for conc exp    
else 
    TrialType=1; %1 for id exp    
end 

%% Get odor repeats
for i = 1:NOdors %for each odor
    reps = find(StimList==StimList(i)); % get indices of all repeats for each odor
    odor_repeats(i) = {reps};
    
end 

%% Classify odor type based on experiment type

switch TrialType
    case 0
disp('concentration exp')
Trial.Id(1,:) = [1 6 11 16]; %Odor 1
Trial.Id(2,:) = [2 7 12 17]; %Odor 2
Trial.Id(3,:) = [3 8 13 18]; %Odor 3
Trial.Id(4,:) = [4 9 14 19]; %Odor 4
Trial.Id(5,:) = [5 10 15 20]; % Odor 5

Trial.Conc(1,:) = [1 2 3 4 5];% concentration 10^-1 
Trial.Conc(2,:) = [6 7 8 9 10];% concentration 10^-2 
Trial.Conc(3,:) = [11 12 13 14 15];% concentration 10^-3 
Trial.Conc(4,:) = [16 17 18 19 20];% concentration 10^-4 



    case 1
disp('identity exp')   
%following lines might not be useful 
% O1 = 1;O2 = 2; O3 = 3;O4 = 4; O5 = 5; O6 = 6; O7 = 7;O8 = 8;
% O9= 9; O10= 10; O11= 11; O12= 12; O13 = 13; O14 = 14; O15 = 15; O16 = 16; 

end 



end