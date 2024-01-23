%FirstvsSecondHalf_ActivevsPassiveReplay_across_animals_script
% Comparing first half and second half of repeats of the active replay to passive replay

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

    [Thismodulated_units,Thismodulation_score] = GetActivePassiveModulatedUnits(SessionName {session}, whichunit, whichARrep, whichPRrep);
    
    Nb_unit(session) = size(Thismodulated_units,1);
    Nb_mod_units_AR(session) = sum(Thismodulated_units(:,1));
    
    perc_modulated_AR(session) = 100*Nb_mod_units_AR(session)/Nb_unit(session);
    
end

%% all units - second halft of active replay repeats

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

   [FiltThismodulated_units,Thismodulation_score] = GetActivePassiveModulatedUnits(SessionName {session}, whichunit, whichARrep, whichPRrep);
    
    FiltNb_unit(session) = size(FiltThismodulated_units,1);
    FiltNb_resp_units_AR(session) = sum(FiltThismodulated_units(:,1));
    
    Filtperc_modulated_AR(session) = 100*FiltNb_resp_units_AR(session)/FiltNb_unit(session);


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

