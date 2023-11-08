function [NumTrials, SuccessRate, MaxTargetHold, TotalTargetStay,...
    ExpertReached, AllMotorLocInTrial, AllCenteredLeverInTrial] ...
    = GetLearningCurveInfo(MyFilePath)

% Extracting the number of trials, success rate, target hold mean time, 
% total target stay time for each session up until the no odor session 
% (= mouse considered expert)
% ExpertReached gives you the total number of session before no odor
% Also getting the motor locations when trial is ON and lever locations
% centered on the tz when trial is ON

if ~exist(MyFilePath,'dir')
    disp('Folder not found')
end 
% session 
cd(MyFilePath)
sessionFiles = dir('*_r0.mat') ; % CHANGE TO ACCOMODATE double sessions a day
numFiles = length(sessionFiles) ; 

% initializing variables
NumTrials = NaN(numFiles,1);
SuccessRate = NaN(numFiles,1);
MaxTargetHold = NaN(numFiles,1);
TotalTargetStay = NaN(numFiles,1);
AllMotorLocInTrial = num2cell(NaN(1,numFiles)); 
AllCenteredLeverInTrial = num2cell(NaN(1,numFiles)); 
SessionCounter =0;


for k = 1 : numFiles
    % Get file name of one behavior mat file
    thisFileName = fullfile(pwd, sessionFiles(k).name);
    fprintf('Processing "%s".\n', sessionFiles(k).name);
    % Load behavior variables
        if  any(strcmp(sessionFiles(k).name,{'S1_20230307_r0.mat',...
                'S6_20230525_r0.mat','S11_20230621_r0.mat'}))
            fprintf('Session skipped because problematic')
            continue % this session is problematic so skip this loop iteration and go to next one
        end
     [TrialInfo,Traces,TargetZones] = PreprocessSmellocatorBehavior(thisFileName);
       
    % Break out of loop if no odor session session 
    if any(strcmp(TrialInfo.Perturbation,'NoOdor'))
        break
    end 
    
    % Getting trials number and sucess rate of this session
    NumTrials(k) = length(TrialInfo.TrialID); 
    SuccessRate(k) = (sum(TrialInfo.Success(:,1))/NumTrials(k))*100;
    
    % Getting Max target hold in tz without exiting
    TrialMaxTargetHold = zeros(NumTrials(k),1);
    for i = 1:NumTrials(k)
        if ~isempty(TrialInfo.InZone{i,1})
            TrialMaxTargetHold(i) = max(TrialInfo.InZone{i,1}(:,2) -TrialInfo.InZone{i,1}(:,1));
        end
    end
    MaxTargetHold(k) = round(mean(TrialMaxTargetHold)*1000);

    % Getting total stay time in tz (with exits and re-entry) 
    TrialTotalTargetStay = zeros(NumTrials(k),1);
    for i = 1:NumTrials(k)
        if ~isempty(TrialInfo.InZone{i,1})
            TrialTotalTargetStay(i) = sum(TrialInfo.InZone{i,1}(:,2) -TrialInfo.InZone{i,1}(:,1));
        end
    end
    TotalTargetStay(k) = round(mean(TrialTotalTargetStay)*1000);

    % Get the distribution of motor locations when trial is on for all trials 
    AllMotorLocations = vertcat(Traces.Motor{:,:});
    AllTrialStates = vertcat(Traces.Trial{:,:});
    AllMotorLocInTrial{1,k} = AllMotorLocations(AllTrialStates~=0);

    %Get the distribution of lever locations centered on tz when trial is
    %on for all trials
    TzPositionThisTrial = NaN(1,NumTrials(k)); %init
    TzCenteredLever = num2cell(NaN(1,NumTrials(k))); %init

    for trial = 1:NumTrials(k)
        %get list of tz position for each trial:
        TzPositionThisTrial(trial) = TargetZones(TrialInfo.TargetZoneType(trial),2);
        %center lever location on tz of that trial:
        TzCenteredLever{1,trial} = Traces.Lever{1,trial}-TzPositionThisTrial(trial);
    end

    AllTzCenteredLever = vertcat(TzCenteredLever{:,:});
    AllTrialStates = vertcat(Traces.Trial{:,:});

    AllCenteredLeverInTrial{1,k} = AllTzCenteredLever(AllTrialStates~=0);

    % Number of sessions analyzed so far
    SessionCounter = SessionCounter+1; 

end
ExpertReached = SessionCounter; % number of session before no odor = expert level
end 

