function [SuccessRateNormal, SuccessRatePerturb, AllCenteredLeverInTrialNormal,AllCenteredLeverInTrialPerturb] =  GetNoOdorSessionInfo(MyFilePath)

[TrialInfo,Traces,TargetZones] = PreprocessSmellocatorBehavior(MyFilePath);

Normaltrials = find(~strcmp(TrialInfo.Perturbation,'NoOdor'));
Perturbtrials = find(strcmp(TrialInfo.Perturbation,'NoOdor'));

%Get Success Rate for both type of trials
SuccessRateNormal = round((sum(TrialInfo.Success(Normaltrials,1))/numel(Normaltrials))*100);
SuccessRatePerturb = round((sum(TrialInfo.Success(Perturbtrials,1))/numel(Perturbtrials))*100);

%Get the distribution of lever locations centered on tz when trial is
%on for normal trials
TzPositionThisTrial = NaN(1,numel(Normaltrials)); %init
TzCenteredLever = num2cell(NaN(1,numel(Normaltrials))); %init

for trial = 1:numel(Normaltrials)
    %get list of tz position for each trial:
    TzPositionThisTrial(trial) = TargetZones(TrialInfo.TargetZoneType(Normaltrials(trial)),2);
    %center lever location on tz of that trial:
    TzCenteredLever{1,trial} = Traces.Lever{1,Normaltrials(trial)}-TzPositionThisTrial(trial);
end

AllTzCenteredLever = vertcat(TzCenteredLever{:,:});
AllTrialStates = vertcat(Traces.Trial{:,Normaltrials(trial)});

AllCenteredLeverInTrialNormal = AllTzCenteredLever(AllTrialStates~=0);

clearvars TzPositionThisTrial TzCenteredLever AllTzCenteredLever AllTrialStates
%Get the distribution of lever locations centered on tz when trial is
%on for perturb trials
TzPositionThisTrial = NaN(1,numel(Perturbtrials)); %init
TzCenteredLever = num2cell(NaN(1,numel(Perturbtrials))); %init

for trial = 1:numel(Perturbtrials)
    %get list of tz position for each trial:
    TzPositionThisTrial(trial) = TargetZones(TrialInfo.TargetZoneType(Perturbtrials(trial)),2);
    %center lever location on tz of that trial:
    TzCenteredLever{1,trial} = Traces.Lever{1,Perturbtrials(trial)}-TzPositionThisTrial(trial);
end

AllTzCenteredLever = vertcat(TzCenteredLever{:,:});
AllTrialStates = vertcat(Traces.Trial{:,Perturbtrials(trial)});

AllCenteredLeverInTrialPerturb = AllTzCenteredLever(AllTrialStates~=0);

end