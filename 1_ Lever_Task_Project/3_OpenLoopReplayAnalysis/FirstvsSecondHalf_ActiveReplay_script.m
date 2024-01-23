%FirstvsSecondHalf_ActiveReplay_across_animals_script
%Comparing 1st half and 2nd half of repeats of the active replay against
%closed loop

% Classic open loop sessions
SessionName = {'O3/O3_20211005_r0_processed.mat',...
'O8/O8_20220702_r0_processed.mat',...
'O9/O9_20220630_r0_processed.mat',...
'S1/S1_20230314_r0_processed.mat',...
'S6/S6_20230727_r0_processed.mat',... 
'S7/S7_20230707_r0_processed.mat',... 
'S12/S12_20230727_r0_processed.mat'...
};
%% 'S3/S3_20230321_r0_processed.mat',... % this mouse is problematic in terms of rec quality so removed 
%% 'S11/S11_20230812_r0_processed.mat',... % not including this mouse because only 4 rep

% O3 10 rep
% O8 10 rep
% O9 10 rep
% S1 6 rep
% S6 7 rep
% S7 8 rep
% S12 8 rep


%% for all units - first half of active replay repeats

ResidualsMean = [];
ResidualsCI95 = [];
ResidualsMedian = [];

whichunit = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    
    if any(strcmp (MouseName{session}, {'O3', 'O8', 'O9'}))
        whichARrep = 1 + (1:4);
        whichPRrep = 1 + 10 + (1:5);
    end
    if strcmp(MouseName{session}, 'S1')
        whichARrep = 1 + (1:3);
        whichPRrep = 1 + 6 + (1:4);
    end
    if strcmp(MouseName{session}, 'S6')
        whichARrep = 1 + (1:3);
        whichPRrep = 1 + 7 + (1:5);
    end
    if strcmp(MouseName{session}, 'S7')
        whichARrep = 1 + (1:4);
        whichPRrep = 1 + 8 + (1:5);
    end
    if strcmp(MouseName{session}, 'S12')
        whichARrep = 1 + (1:4);
        whichPRrep = 1 + 8 + (1:4);
    end

    [Thismodulated_units,Thismodulation_score, ThisResidualsMean,...
        ThisResidualsMedian, ThisResidualsCI95] = GetReplayModulatedUnits_v2(SessionName {session}, whichunit, whichARrep, whichPRrep);
    
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

%% all units - second halft of active replay repeats

FiltResidualsMean = [];
FiltResidualsCI95 = [];
FiltResidualsMedian = [];

whichunit = [];

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    
    if any(strcmp (MouseName{session}, {'O3', 'O8', 'O9'}))
        whichARrep = 1 + (7:10);
        whichPRrep = 1 + numel(whichARrep) + (1:5);
    end
    if strcmp(MouseName{session}, 'S1')
        whichARrep = 1 + (4:6);
        whichPRrep = 1 + 6 + (1:4);
    end
    if strcmp(MouseName{session}, 'S6')
        whichARrep = 1 + (5:7);
        whichPRrep = 1 + 7 + (1:5);
    end
    if strcmp(MouseName{session}, 'S7')
        whichARrep = 1 + (5:8);
        whichPRrep = 1 + 8 + (1:5);
    end
    if strcmp(MouseName{session}, 'S12')
        whichARrep = 1 + (5:8);
        whichPRrep = 1 + numel(whichARrep) + (1:4);
    end 

    [FiltThismodulated_units, FiltThismodulation_score, FiltThisResidualsMean,...
        FiltThisResidualsMedian, FiltThisResidualsCI95] = GetReplayModulatedUnits_v2(SessionName {session}, whichunit, whichARrep, whichPRrep);
    
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


%%  PLOTS for all units for middle PR %%
%% percentage modulated by AR percentage modulated by middle PR
figure(1) 
subplot(1,2,1) % for each mouse
y = [perc_modulated_AR; Filtperc_modulated_AR]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];
ylim([0 100])
yticks(0:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('1st half active replay', '2nd half active replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_AR ; Filtperc_modulated_AR]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'1st half active replay', '2nd half active replay'})

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