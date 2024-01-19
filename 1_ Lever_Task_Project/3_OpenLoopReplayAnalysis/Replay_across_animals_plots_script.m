%Replay_across_animals_plots_script


% Specifify which type of OL session
%OL = "classic" ;
OL = "free";

switch OL
    case "classic"
% Classic open loop sessions
SessionName = {'O3/O3_20211005_r0_processed.mat',...
'O8/O8_20220702_r0_processed.mat',...
'O9/O9_20220630_r0_processed.mat',...
'S1/S1_20230314_r0_processed.mat',...
'S6/S6_20230727_r0_processed.mat',... 
'S7/S7_20230707_r0_processed.mat',... 
'S11/S11_20230812_r0_processed.mat',...
'S12/S12_20230727_r0_processed.mat'};
%% 'S3/S3_20230321_r0_processed.mat',... % this mouse is problematic in terms of rec quality so removed 

    case "free"
% Free lever Passive Replay
SessionName = {...  
'S1/S1_20230403_r0_processed.mat',... 
'S6/S6_20230718_r0_processed.mat',...
'S7/S7_20230622_r0_processed.mat',... 
'S11/S11_20230805_r0_processed.mat',...
'S12/S12_20230804_r0_processed.mat'...
};

end
%% for all units

ResidualsMean = [];
ResidualsCI95 = [];
ResidualsMedian = [];

whichunit = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    [Thismodulated_units,Thismodulation_score, ThisResidualsMean, ...
        ThisResidualsMedian, ThisResidualsCI95] = GetReplayModulatedUnits_v2(SessionName {session}, whichunit);
    
    Nb_unit(session) = size(Thismodulated_units,1);
    Nb_mod_units_AR(session) = sum(Thismodulated_units(:,1));
    Nb_mod_units_PR(session) = sum(Thismodulated_units(:,2));
    
    perc_modulated_AR(session) = 100*Nb_mod_units_AR(session)/Nb_unit(session);
    perc_modulated_PR(session) = 100*Nb_mod_units_PR(session)/Nb_unit(session);
    
    if session ==1
        ResidualsMean = ThisResidualsMean;
        ResidualsCI95 = ThisResidualsCI95;
        ResidualsMedian = ThisResidualsMedian;
    elseif session >1
        ResidualsMean = cat(1,ResidualsMean, ThisResidualsMean);
        ResidualsCI95 = cat(1,ResidualsCI95, ThisResidualsCI95);
        ResidualsMedian = cat(1,ResidualsMedian, ThisResidualsMedian);
    end

end

%% filtered by responsive units

FiltResidualsMean = [];
FiltResidualsCI95 = [];
FiltResidualsMedian = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});

    [responsive_units, ~] = GetResponsiveUnits(SessionName {session});
    whichunit = find(responsive_units~=0);

    [FiltThismodulated_units, FiltThismodulation_score, FiltThisResidualsMean, ...
        FiltThisResidualsMedian, FiltThisResidualsCI95] = GetReplayModulatedUnits_v2(SessionName {session}, whichunit);
    
    FiltNb_unit(session) = size(FiltThismodulated_units,1);
    FiltNb_resp_units_AR(session) = sum(FiltThismodulated_units(:,1));
    FiltNb_resp_units_PR(session) = sum(FiltThismodulated_units(:,2));
    
    Filtperc_modulated_AR(session) = 100*FiltNb_resp_units_AR(session)/FiltNb_unit(session);
    Filtperc_modulated_PR(session) = 100*FiltNb_resp_units_PR(session)/FiltNb_unit(session);

    if session ==1
        FiltResidualsMean = FiltThisResidualsMean;
        FiltResidualsCI95 = FiltThisResidualsCI95;
        FiltResidualsMedian = FiltThisResidualsMedian;
    elseif session>1
        FiltResidualsMean = cat(1,FiltResidualsMean, FiltThisResidualsMean);
        FiltResidualsCI95 = cat(1,FiltResidualsCI95, FiltThisResidualsCI95);
        FiltResidualsMedian = cat(1,FiltResidualsMedian, FiltThisResidualsMedian);
    end

end


%%  PLOTS for all units %%
%% percentage modulated by AR percentage modulated by PR
figure(1) 
subplot(1,2,1) % for each mouse
y = [perc_modulated_AR; perc_modulated_PR]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];
ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('active replay', 'passive replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_AR ; perc_modulated_PR]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'active replay','passive replay'})

%% figure with just passive replay
figure(2) 
subplot(1,2,1) % for each mouse
y = [perc_modulated_PR]';
b = bar(y);
b.FaceColor = [0 0 0];
ylim([0 100])
yticks(0:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')

legend('passive replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_PR]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'passive replay'})

%% scatter plot of self vs. across condition residuals
 %[ResidualsMean{1,1};ResidualsMean{2,1}]
figure(3) 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
    for i = 1:size(ResidualsMedian,1)
        MedianResiduals = [MedianResiduals; ResidualsMedian{i,whichtype}];
        STDResiduals = [STDResiduals; ResidualsCI95{i,whichtype}];
    end
end

    subplot(1,2,1) %CL-PR with respect to PR-PR
    [~,sortorder] = sort(MedianResiduals(:,5));
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) + STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) - STDResiduals(sortorder,5),':k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5),'k');
    plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,2),'or');
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);
    
    subplot(1,2,2) %CL-OL with respect to OL-OL
    [~,sortorder] = sort(MedianResiduals(:,3)+STDResiduals(:,3));
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) + STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) - STDResiduals(sortorder,3),':k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3),'k');
    plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,1),'or');
    axis square;
    set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);


%% PLOTS for filtered units %%
%% percentage modulated by AR percentage modulated by PR
figure(4) 
subplot(1,2,1) % for each mouse
y = [Filtperc_modulated_AR; Filtperc_modulated_PR]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];
ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('active replay', 'passive replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [Filtperc_modulated_AR ; Filtperc_modulated_PR]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'active replay', 'passive replay'})


%% scatter plot of self vs. across condition residuals for the filtered data
figure(5) 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
    for i = 1:size(FiltResidualsMedian,1)
        MedianResiduals = [MedianResiduals; FiltResidualsMedian{i,whichtype}];
        STDResiduals = [STDResiduals; FiltResidualsCI95{i,whichtype}];
    end
end

subplot(1,2,1) %CL-PR with respect to PR-PR
[~,sortorder] = sort(MedianResiduals(:,5));
plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) + STDResiduals(sortorder,5),':k');
plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5) - STDResiduals(sortorder,5),':k');
x = MedianResiduals(sortorder,5);
y1 = MedianResiduals(sortorder,5) + STDResiduals(sortorder,5);
y2 = MedianResiduals(sortorder,5) - STDResiduals(sortorder,5);
patch(y1, y2, 'g')
plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,5),'k');
plot(MedianResiduals(sortorder,5),MedianResiduals(sortorder,2),'or');
axis square;
set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);

subplot(1,2,2) %CL-OL with respect to OL-OL
[~,sortorder] = sort(MedianResiduals(:,3)+STDResiduals(:,3));
plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) + STDResiduals(sortorder,3),':k');
plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3) - STDResiduals(sortorder,3),':k');
plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,3),'k');
plot(MedianResiduals(sortorder,3),MedianResiduals(sortorder,1),'or');
axis square;
set(gca,'TickDir','out','YLim',[0 20],'XLim',[0 20]);


% %% TESTING PLOTS FOR DISTRIBUTION
% %% distribution of mean
% MedianResiduals =[]; STDResiduals =[];
% for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
%     for i = 1:size(ResidualsMean,1)
%         MedianResiduals = [MedianResiduals; ResidualsMean{i,whichtype}];
%         STDResiduals = [STDResiduals; ResidualsCI95{i,whichtype}];
%     end
% end
% figure(5)
% [N,edges] = histcounts(MedianResiduals(sortorder,5));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(sortorder,2));
% plot(N/sum(N))
% 
% figure(6)
% [N,edges] = histcounts(MedianResiduals(sortorder,3));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(sortorder,1));
% plot(N/sum(N))
% 
% %% distribution of median
% MedianResiduals =[]; STDResiduals =[];
% for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
%     for i = 1:size(ResidualsMean,1)
%         MedianResiduals = [MedianResiduals; ResidualsMedian{i,whichtype}];
%         STDResiduals = [STDResiduals; ResidualsCI95{i,whichtype}];
%     end
% end
% 
% figure(7)
% [N,edges] = histcounts(MedianResiduals(sortorder,5));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(sortorder,2));
% plot(N/sum(N))
% 
% figure(8)
% [N,edges] = histcounts(MedianResiduals(sortorder,3));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(sortorder,1));
% plot(N/sum(N))
% 
% %% for odor 1
% MedianResiduals =[]; STDResiduals =[];
% whichtype = 2; % = PSTHs comparisons for odor 1 to 3
%     for i = 1:size(ResidualsMean,1)
%         MedianResiduals = [MedianResiduals; ResidualsMedian{i,whichtype}];
%         STDResiduals = [STDResiduals; ResidualsCI95{i,whichtype}];
%     end
% 
% 
% figure(8)
% [N,edges] = histcounts(MedianResiduals(:,5));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(:,2));
% plot(N/sum(N))
% 
% figure(9)
% [N,edges] = histcounts(MedianResiduals(:,3));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(:,1));
% plot(N/sum(N))
% 
% %% for one mouse
% 
% MedianResiduals =[]; STDResiduals =[];
% whichtype = 2; % = PSTHs comparisons for odor 1 to 3
%     for i = 1:size(ThisResidualsMean,1)
%         MedianResiduals = [MedianResiduals; ThisResidualsMedian{i,whichtype}];
%         STDResiduals = [STDResiduals; ThisResidualsCI95{i,whichtype}];
%     end
% 
% 
% figure(10)
% [N,edges] = histcounts(MedianResiduals(:,5));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(:,2));
% plot(N/sum(N))
% 
% figure(11)
% [N,edges] = histcounts(MedianResiduals(:,3));
% plot(N/sum(N))
% hold on
% [N,edges] = histcounts(MedianResiduals(:,1));
% plot(N/sum(N))