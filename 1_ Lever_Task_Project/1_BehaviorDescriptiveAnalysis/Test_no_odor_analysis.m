%Test_no_odor_analysis

Session = '/Users/mariedussauze/Desktop/Analysis/data/Behavior/S1/S1_20230312_r0.mat';
[TrialInfo,Traces,TargetZones] = PreprocessSmellocatorBehavior(Session);

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

%% PLOTTING %%

fig1 = figure(1); % histograms of lever location in trials for all sessions
BinWidth = 0.6;
bins = -4.5:BinWidth:4.5;% Bin width of .6 allows to have entire tz for lever as one of the bins
%TZCenterInd = find(BinCenters==0);
BinCenters = (bins(1:end-1)+bins(2:end))/2;
fwhmx = NaN(1,1);
[N,edges] = histcounts(AllCenteredLeverInTrialNormal, bins);
%h = histogram('BinEdges',edges,'BinCounts',(N/sum(N))*100);
occupancy = N/sum(N);
plot(BinCenters,(occupancy*100))


fig2 = figure(2); % histograms of lever location in trials for all sessions
BinWidth = 0.6;
bins = -4.5:BinWidth:4.5;% Bin width of .6 allows to have entire tz for lever as one of the bins
%TZCenterInd = find(BinCenters==0);
BinCenters = (bins(1:end-1)+bins(2:end))/2;
fwhmx = NaN(1,1);
[N,edges] = histcounts(AllCenteredLeverInTrialPerturb, bins);
%h = histogram('BinEdges',edges,'BinCounts',(N/sum(N))*100);
occupancy = N/sum(N);
plot(BinCenters,(occupancy*100))
area(BinCenters,(occupancy*100), 'FaceColor', [0.5 0.5 0.5])




