function [SuccessRateNormal, SuccessRatePerturb] =  GetOffsetSessionInfo(MyFilePath)

[TrialInfo,Traces,TargetZones] = PreprocessSmellocatorBehavior(MyFilePath);

Normaltrials = find(~strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template'));
Perturbtrials = find(strcmp(TrialInfo.Perturbation(:,1),'Offset-II-Template'));

%Get Success Rate for both type of trials
SuccessRateNormal = round((sum(TrialInfo.Success(Normaltrials,1))/numel(Normaltrials))*100);
SuccessRatePerturb = round((sum(TrialInfo.Success(Perturbtrials,1))/numel(Perturbtrials))*100);

end

