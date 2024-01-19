%No_odor_across_mice_script

SessionName = {'S1/S1_20230312_r0.mat', 'S3/S3_20230316_r0.mat', 'S6/S6_20230622_r0.mat',...
    'S7/S7_20230530_r0.mat', 'S11/S11_20230724_r0.mat', 'S12/S12_20230724_r0.mat'};
BehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Behavior';

CenteredLeverInTrialNormal_all = [];
CenteredLeverInTrialPerturb_all = [];

for session = 1:length(SessionName)
    MyFilePath = fullfile(BehaviorPath,SessionName{session});
[SuccessRateNormal, SuccessRatePerturb, AllCenteredLeverInTrialNormal,...
    AllCenteredLeverInTrialPerturb] =  GetNoOdorSessionInfo(MyFilePath);

SuccessRateNormal_all(session) = SuccessRateNormal;
SuccessRatePerturb_all(session) = SuccessRatePerturb;

if strcmp(SessionName{session}, 'S1/S1_20230312_r0.mat')
CenteredLeverInTrialNormal_all = vertcat(CenteredLeverInTrialNormal_all,AllCenteredLeverInTrialNormal) ;
CenteredLeverInTrialPerturb_all = vertcat(CenteredLeverInTrialPerturb_all, AllCenteredLeverInTrialPerturb);
end

end

%% PLOTTING - TRYING TO CHARACTERIZE THE SPREAD OF LEVER DISTRIBUTION FOR BOTH CONDITIONS
BinWidth = 0.45;
bins = -4.5:BinWidth:4.5;% Bin width of .6 allows to have entire tz for lever as one of the bins
%TZCenterInd = find(BinCenters==0);
BinCenters = (bins(1:end-1)+bins(2:end))/2;

fig1 = figure(1); % histograms of lever location in trials for all sessions
[N,edges] = histcounts(CenteredLeverInTrialNormal_all, bins);
%h = histogram('BinEdges',edges,'BinCounts',(N/sum(N))*100);
occupancy_normal = N/sum(N);
plot(BinCenters,(occupancy_normal*100))
area(BinCenters,(occupancy_normal*100))

fig2 = figure(2); % histograms of lever location in trials for all sessions
[N,edges] = histcounts(CenteredLeverInTrialPerturb_all, bins);
%h = histogram('BinEdges',edges,'BinCounts',(N/sum(N))*100);
occupancy_perturb = N/sum(N);
plot(BinCenters,(occupancy_perturb*100))
area(BinCenters,(occupancy_perturb*100), 'FaceColor', [0.7 0.7 0.7])

figure(3)
area(BinCenters,(occupancy_normal*100))
hold on
area(BinCenters,(occupancy_perturb*100), 'FaceColor', [0.7 0.7 0.7])

figure(4)
plot(BinCenters,(occupancy_normal*100), 'LineWidth',2)
hold on
plot(BinCenters,(occupancy_perturb*100), 'LineWidth',2)
xline(0)

%% PLOTS
figure(5) % Success Rate
for i = 1:length(SessionName)
    plot([1 2], [SuccessRateNormal_all(i) SuccessRatePerturb_all(i)], '.-',...
        'LineWidth',2, 'MarkerSize',12)
    set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
    xlim([0 3])
    ylim([0 100])
    xticks([1 2])
    yticks(20:20:100)
    xticklabels({'Normal', 'No Odor'})
    ylabel('Success Rate (%)')
    hold on
end
