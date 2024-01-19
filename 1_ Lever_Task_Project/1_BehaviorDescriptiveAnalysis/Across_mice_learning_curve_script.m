%Across_mice_learning_curve_script

MouseName = {'S1', 'S3', 'S6', 'S7', 'S11', 'S12'};
BehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Behavior';

for mouse = 1:length(MouseName)
    MyFilePath = fullfile(BehaviorPath,MouseName{mouse});
    [NumTrials, SuccessRate, MaxTargetHold, TotalTargetStay, TrialLength, ...
    ExpertReached, AllMotorLocInTrial, AllCenteredLeverInTrial]...
    = GetLearningCurveInfo(MyFilePath);

    NumTrials_expert(mouse,:) = NumTrials(ExpertReached-2:ExpertReached); % keeping only the last 3 baseline sessions before no odor
    SuccessRate_expert(mouse,:) = SuccessRate(ExpertReached-2:ExpertReached);
    TotalTargetStay_expert(mouse,:) = TotalTargetStay(ExpertReached-2:ExpertReached);
    TrialLength_expert(mouse,:) = TrialLength(ExpertReached-2:ExpertReached);
    MaxTargetHold_expert(mouse,:) = MaxTargetHold(ExpertReached-2:ExpertReached);

    ExpertReached_all(mouse,:) = ExpertReached;
    NumTrials_all{mouse} = NumTrials; 
    SuccessRate_all{mouse} = SuccessRate;
    MaxTargetHold_all{mouse} = MaxTargetHold;
    TotalTargetStay_all{mouse} = TotalTargetStay;
    TrialLength_all{mouse} = TrialLength;
end


%% CORRECTING FOR MISSING INFO BASED ON LOGS
% S1 20230307
MaxTargetHold_all{1}(19) = 335;
SuccessRate_all{1}(19) = 80;
% S6 20230525
MaxTargetHold_all{3}(8) = 162;
SuccessRate_all{3}(8) = 63;
% S11 20230621
MaxTargetHold_all{5}(7) = 259.5;
SuccessRate_all{5}(7) = 81;

%% PLOTTING %%
%% Target hold plots
figure(1) % Leaning curve of max target hold 
subplot(1,4,[1 2])
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

subplot(1,4,3) %how many sessions to be expert 
y = reshape(MaxTargetHold_expert,numel(MaxTargetHold_expert),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 500])
ylabel('Max target hold (ms)')
yticks(100:100:500)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')

subplot(1,4,4) %how many sessions to be expert 
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

%% Number of trials
figure(2) %num trials across mice and their 3 last stabilizing sessions
subplot(1,3,[1 2]) %num trials across sessions
for i=1:length(MouseName)
plot(NumTrials_all{i}, '.-','LineWidth',2, 'MarkerSize',12)
hold on
ylim([0 500])
xlabel('Session Number')
ylabel('Number of trials')
yticks(100:100:500)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
end

subplot(1,3,3) %num trials across mice and their 3 last stabilizing sessions
y = reshape(NumTrials_expert,numel(NumTrials_expert),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 500])
ylabel('Nb of trials per session')
yticks(0:100:500)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')

%% Success rate 
figure(3) 
subplot(1,3,[1 2]) %success rate across sessions
for i=1:length(MouseName)
plot(SuccessRate_all{i}, '.-','LineWidth',2, 'MarkerSize',12)
hold on
ylim([0 100])
xlabel('Session Number')
ylabel('Success rate')
yticks([0 50 100])
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
end

subplot(1,3,3) %success rate across mice and their 3 last stabilizing sessions
y = reshape(SuccessRate_expert,numel(SuccessRate_expert),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
ylim([0 100])
ylabel('Success rate (%')
yticks([0 50 100])
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')

%% Total Target Stay divided by trial length - get a sense of how much time they spend in TZ vs the rest
figure(4) % Leaning curve of tot stay in tz
subplot(1,3,[1 2])
for i=1:length(MouseName)
plot(TotalTargetStay_all{i}./TrialLength_all{i}, '.-','LineWidth',2, 'MarkerSize',12)
hold on
ylim([0 0.5])
xlabel('Session Number')
ylabel('Max target hold (ms)/Trial length')
yticks([0 0.5])
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
end

%% Length of trial
figure(5) % Leaning curve of tot stay in tz
subplot(1,3,[1 2])
subplot(1,3,3) %success rate across mice and their 3 last stabilizing sessions
y = reshape(TrialLength_expert,numel(TrialLength_expert),1);
boxchart(y,'BoxFaceColor',"#A2142F", 'MarkerColor',"#A2142F", 'BoxWidth',0.4)
xticks([])
%ylim([0 100])
ylabel('Trial Length (ms)')
%yticks([0 100])
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')