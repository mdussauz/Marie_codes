%Offset_Across_Mice_Success_script

SessionName = {'S1/S1_20230405_r0.mat', 'S3/S3_20230405_r0.mat', 'S6/S6_20230721_r0.mat',...
    'S7/S7_20230616_r0.mat', 'S11/S11_20230807_r0.mat', 'S12/S12_20230806_r0.mat'};
BehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Behavior';

CenteredLeverInTrialNormal_all = [];
CenteredLeverInTrialPerturb_all = [];

for session = 1:length(SessionName)
    MyFilePath = fullfile(BehaviorPath,SessionName{session});
    [SuccessRateNormal, SuccessRatePerturb] =  GetOffsetSessionInfo(MyFilePath);

    SuccessRateNormal_all(session) = SuccessRateNormal;
    SuccessRatePerturb_all(session) = SuccessRatePerturb;
end

%% PLOTS %%
%%
figure % Success Rate
for i = 1:length(SessionName)
    plot([1 2], [SuccessRateNormal_all(i) SuccessRatePerturb_all(i)], '.-',...
        'LineWidth',2, 'MarkerSize',12)
    set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
    xlim([0 3])
    ylim([0 100])
    xticks([1 2])
    yticks(25:25:100)
    xticklabels({'Normal trials', 'Offset trials'})
    ylabel('Success Rate (%)')
    hold on
end