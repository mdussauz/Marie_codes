% HaltFlip_across_animals_script

%% Decide whether to filter by responsiveness or not
filter_by_responsiveness =0;

%% Extract data
% AON
SessionName = {'O8/O8_20220704_r0_processed.mat',...
    'O9/O9_20220702_r1_processed.mat',...
    'S1/S1_20230327_r0_processed.mat',...
    'S6/S6_20230710_r0_processed.mat', ...
    'S7/S7_20230608_r0_processed.mat',...
    'S11/S11_20230801_r0_processed.mat',...
    'S12/S12_20230731_r0_processed.mat'};

% APC
% SessionName = {'Q4/Q4_20221109_r0_processed.mat',...
%     'Q8/Q8_20221207_r0_processed.mat',...
%     'Q9/Q9_20221116_r0_processed.mat'};



for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    if filter_by_responsiveness
        [responsive_units, ~] = GetResponsiveUnits(SessionName {session});
        whichunit = find(responsive_units~=0);
        [mirror_modulation, location_modulation, replay_modulation, tuning_modulation] = GetHaltModulatedUnits(SessionName{session},whichunit);
    else
        [mirror_modulation, location_modulation, replay_modulation, tuning_modulation] = GetHaltModulatedUnits(SessionName{session});
    end
    Nb_unit(session) = size(location_modulation,2);
    Nb_mod_units_loc(session) = sum(location_modulation~=0);
    Nb_mod_units_replay(session) = sum(replay_modulation~=0);
    Nb_mod_units_tuning(session) = sum(tuning_modulation~=0);
    
    % Cells that change across both location tuning and replay
    overlap_loc_replay(session) = numel(intersect(find(location_modulation~=0), find(replay_modulation~=0)));
    % Cells that change across both passive tuning and replay
    overlap_tuning_replay(session) = numel(intersect(find(tuning_modulation~=0), find(replay_modulation~=0)));

    % Excited/inhibited 
    exc_loc = 1*(location_modulation==1); %1 if excited by any event 
    inh_loc = 2*(location_modulation==2); %1 if inihibited by any event 
    perc_exc_loc(session) = 100*sum(exc_loc)/Nb_unit(session);
    perc_inh_loc(session) = 100*sum(inh_loc)/Nb_unit(session);

    exc_tun = 1*(tuning_modulation==1); %1 if excited by any event 
    inh_tun = 2*(tuning_modulation==2); %1 if inihibited by any event 
    perc_exc_tun(session) = 100*sum(exc_tun)/Nb_unit(session);
    perc_inh_tun(session) = 100*sum(inh_tun)/Nb_unit(session);

    
    %Calculate percentages
    perc_modulated_loc(session) = 100*Nb_mod_units_loc(session)/Nb_unit(session);
    perc_modulated_replay(session) = 100*Nb_mod_units_replay(session)/Nb_unit(session);
    perc_modulated_tuning(session) = 100*Nb_mod_units_tuning(session)/Nb_unit(session);

    perc_overlap_loc(session) = 100*overlap_loc_replay(session)/Nb_unit(session);
    perc_overlap_tun(session) = 100*overlap_tuning_replay(session)/Nb_unit(session);

end 


%%  PLOTS for all units %%
%% percentage modulated compared to location tuning curve, replay, passive tuning
figure(1) 
subplot(1,2,1) % for each mouse
y = [perc_modulated_loc; perc_modulated_replay; perc_modulated_tuning]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

legend('tuning curve', 'halt replay', 'passive tuning')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)

subplot(1,2,2) % across mice
y = [perc_modulated_loc; perc_modulated_replay; perc_modulated_tuning]';
boxchart(y)
hold on
plot(mean(y), '.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'tuning curve', 'halt replay', 'passive tuning'})
legend('data','mean')
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)

%% Putative sensorimotor Overlap of replay halt and location modulated
% as well as overlap of replay and passive tuning

figure(2)
subplot(1,2,1)
y = [perc_overlap_loc; perc_overlap_tun]'; 
b = bar(y);
b.FaceColor = [0.5 0.5 0.5];
ylim([0 100])
yticks(0:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2)
boxchart(y)
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'putative loc sensorimotor', 'putative tunung sensorimotor'})

%% tuning only
figure(3)
subplot(1,2,1) % for each mouse
y = [perc_modulated_loc]';
b = bar(y);
b.FaceColor = [0.5 0.5 0.5];

ylim([0 100])
yticks(0:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('tuning curve')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_loc]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y), '.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'tuning curve'})

%% passive halt only
figure(4)
subplot(1,2,1) % for each mouse
y = [perc_modulated_replay]';
b = bar(y);
b.FaceColor = [0.5 0.5 0.5];

ylim([0 100])
yticks(0:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)

legend('passive halt')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_replay]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(0:20:100)
hold on
plot(mean(y), '.','MarkerSize',12,'Color','k')
hold off
ylabel ('% modulated units')
xticklabels({'passive halt'})

%% excited inhibited compared to location tuning
figure(5) % percentage responsive: excited + inhibited + mix
subplot(1,2,1) % for each mouse
y = [perc_exc_loc; perc_inh_loc]';
b = bar(y,'stacked');
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0.5 0.5];


ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
legend('excited', 'inhibited')
ylabel ('% responsive units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_exc_loc ; perc_inh_loc]';
boxchart(y,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% responsive units')
xticklabels({'excited','inhibited'})

%% excited inhibited compared to passive tuning
figure(6) % percentage responsive: excited + inhibited + mix
subplot(1,2,1) % for each mouse
y = [perc_exc_loc; perc_inh_loc]';
b = bar(y,'stacked');
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [0.5 0.5 0.5];


ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
legend('excited', 'inhibited')
ylabel ('% responsive units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_exc_loc ; perc_inh_loc]';
boxchart(y,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% responsive units')
xticklabels({'excited','inhibited'})
