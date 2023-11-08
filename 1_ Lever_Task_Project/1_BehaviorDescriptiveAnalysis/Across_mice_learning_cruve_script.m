%Across_mice_learning_curve_script

MouseName = {'S1', 'S3', 'S6', 'S7', 'S11', 'S12'};
BehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Behavior';

for mouse = 1:length(MouseName)
    MyFilePath = fullfile(BehaviorPath,MouseName{mouse});
    [NumTrials, SuccessRate, MaxTargetHold, TotalTargetStay,  ExpertReached,...
    AllMotorLocInTrial, AllCenteredLeverInTrial]...
    = GetLearningCurveInfo(MyFilePath);

    NumTrials_all(mouse,:) = NumTrials(ExpertReached-2:ExpertReached); % keeping only the last 3 baseline sessions before no odor
    SuccessRate_all(mouse,:) = SuccessRate(ExpertReached-2:ExpertReached);
    ExpertReached_all(mouse,:) = ExpertReached;
    MaxTargetHold_all{mouse} = MaxTargetHold;
end


%% CORRECTING FOR MISSING INFO BASED ON LOGS
%S1 20230307
 MaxTargetHold_all{1}(19) = 335;
 % S6 20230525
MaxTargetHold_all{3}(8) = 162;
 %S11 20230621
 MaxTargetHold_all{5}(7) = 259.5;
%% PLOTTING %%
figure(1) % Leaning curve of max target hold 
for i=1:length(MouseName)
plot(MaxTargetHold_all{i}, '.-','LineWidth',2, 'MarkerSize',12)
hold on
ylim([0 500])
xlabel('Session Number')
ylabel('Max target hold (ms)')
yticks(100:100:500)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
end

figure(2) %how many sessions to be expert 
y = ExpertReached_all;
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 30])
ylabel('Nb of sessions to expert')
yticks(0:10:30)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')

figure(3) %num trials across mice and their 3 last stabilizing sessions
y = reshape(NumTrials_all,numel(NumTrials_all),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 500])
ylabel('Nb of trials per session')
yticks(0:100:500)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')

figure(4) %success rate across mice and their 3 last stabilizing sessions
y = reshape(SuccessRate_all,numel(SuccessRate_all),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 100])
ylabel('Success rate (%')
yticks([0 100])
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
