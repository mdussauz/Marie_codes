%Replay_across_animals_plots_script
%% Classic open loop sessions
% SessionName = {'O3/O3_20211005_r0_processed.mat',...
% 'O8/O8_20220702_r0_processed.mat',...
% 'O9/O9_20220630_r0_processed.mat',...
% 'S1/S1_20230314_r0_processed.mat',...
% 'S3/S3_20230321_r0_processed.mat',...
% 'S6/S6_20230718_r0_processed.mat',... % !!!this session is a free lever one - to be changed when classic is fixed!!!
% 'S7/S7_20230707_r0_processed.mat',... 
% 'S11/S11_20230812_r0_processed.mat',...
% 'S12/S12_20230727_r0_processed.mat'};

%% Free lever Passive Replay
SessionName = {...   
'S6/S6_20230718_r0_processed.mat',...
'S7/S7_20230622_r0_processed.mat',... 
'S11/S11_20230805_r0_processed.mat'...
};


%% for all units

ResidualsMean = [];
ResidualsCI95 = [];

whichunit = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    [Thismodulated_units,Thismodulation_score, ThisResidualsMean, ThisResidualsCI95] = GetReplayModulatedUnits(SessionName {session}, whichunit);
    
    Nb_unit(session) = size(Thismodulated_units,1);
    Nb_mod_units_AR(session) = sum(Thismodulated_units(:,1));
    Nb_mod_units_PR(session) = sum(Thismodulated_units(:,2));
    
    perc_modulated_AR(session) = 100*Nb_mod_units_AR(session)/Nb_unit(session);
    perc_modulated_PR(session) = 100*Nb_mod_units_PR(session)/Nb_unit(session);
    
    if session ==1
        ResidualsMean = ThisResidualsMean;
        ResidualsCI95 = ThisResidualsCI95;
    elseif session >1
        ResidualsMean = cat(1,ResidualsMean, ThisResidualsMean);
        ResidualsCI95 = cat(1,ResidualsCI95, ThisResidualsCI95);
    end

end

%% filtered by responsive units

FiltResidualsMean = [];
FiltResidualsCI95 = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});

    [responsive_units, ~] = GetResponsiveUnits(SessionName {session});
    whichunit = find(responsive_units~=0);

    [FiltThismodulated_units, FiltThismodulation_score, FiltThisResidualsMean, FiltThisResidualsCI95] = GetReplayModulatedUnits(SessionName {session}, whichunit);
    
    FiltNb_unit(session) = size(FiltThismodulated_units,1);
    FiltNb_resp_units_AR(session) = sum(FiltThismodulated_units(:,1));
    FiltNb_resp_units_PR(session) = sum(FiltThismodulated_units(:,2));
    
    Filtperc_modulated_AR(session) = 100*FiltNb_resp_units_AR(session)/FiltNb_unit(session);
    Filtperc_modulated_PR(session) = 100*FiltNb_resp_units_PR(session)/FiltNb_unit(session);

    if session ==1
        FiltResidualsMean = FiltThisResidualsMean;
        FiltResidualsCI95 = FiltResidualsCI95;
    elseif session>1
        FiltResidualsMean = cat(1,FiltResidualsMean, FiltThisResidualsMean);
        FiltResidualsCI95 = cat(1,FiltResidualsCI95, FiltThisResidualsCI95);
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

legend('active replay', 'passive replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_AR ; perc_modulated_PR]';
boxchart(y)
hold on
plot(mean(y), 'o')
hold off
ylabel ('% modulated units')
xticklabels({'active replay','passive replay'})
legend('data','mean')

%% scatter plot of self vs. across condition residuals
 %[ResidualsMean{1,1};ResidualsMean{2,1}]
figure(2) 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
    for i = 1:size(ResidualsMean,1)
        MedianResiduals = [MedianResiduals; ResidualsMean{i,whichtype}];
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
figure(3) 
subplot(1,2,1) % for each mouse
y = [Filtperc_modulated_AR; Filtperc_modulated_PR]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

legend('active replay', 'passive replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [Filtperc_modulated_AR ; Filtperc_modulated_PR]';
boxchart(y)
hold on
plot(mean(y), 'o')
hold off
ylabel ('% modulated units')
xticklabels({'active replay', 'passive replay'})
legend('data','mean')

%% scatter plot of self vs. across condition residuals
figure(4) 
subplot(1,2,1), hold on
subplot(1,2,2), hold on
MedianResiduals =[]; STDResiduals =[];
for whichtype = 2:4 % = PSTHs comparisons for odor 1 to 3
    for i = 1:size(FiltResidualsMean,1)
        MedianResiduals = [MedianResiduals; FiltResidualsMean{i,whichtype}];
        STDResiduals = [STDResiduals; FiltResidualsCI95{i,whichtype}];
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