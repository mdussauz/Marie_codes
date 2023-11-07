% HaltFlip_across_animals_script

%% Extract data
SessionName = {'O8/O8_20220704_r0_processed.mat',...
    'O9/O9_20220702_r1_processed.mat',...
    'S1/S1_20230327_r0_processed.mat',...
    'S6/S6_20230710_r0_processed.mat', ...
    'S7/S7_20230608_r0_processed.mat',...
    'S11/S11_20230801_r0_processed.mat',...
    'S12/S12_20230731_r0_processed.mat'};


for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});

    [mirror_modulation, location_modulation, replay_modulation, tuning_modulation] = GetHaltModulatedUnits(SessionName{session});

    Nb_unit(session) = size(location_modulation,2);
    Nb_mod_units_loc(session) = sum(location_modulation~=0);
    Nb_mod_units_replay(session) = sum(replay_modulation~=0);
    Nb_mod_units_tuning(session) = sum(tuning_modulation~=0);

    overlap_loc_replay(session) = numel(intersect(find(location_modulation~=0), find(replay_modulation~=0)));
    
    perc_modulated_loc(session) = 100*Nb_mod_units_loc(session)/Nb_unit(session);
    perc_modulated_replay(session) = 100*Nb_mod_units_replay(session)/Nb_unit(session);
    perc_modulated_tuning(session) = 100*Nb_mod_units_tuning(session)/Nb_unit(session);

    perc_overlap(session) = 100*overlap_loc_replay(session)/Nb_unit(session);

end 


%%  PLOTS for all units %%
%% percentage modulated compared to location tuning curve, replay, tuning
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

subplot(1,2,2) % across mice
y = [perc_modulated_loc; perc_modulated_replay; perc_modulated_tuning]';
boxchart(y)
hold on
plot(mean(y), 'o')
hold off
ylabel ('% modulated units')
xticklabels({'tuning curve', 'halt replay', 'passive tuning'})
legend('data','mean')


%% same but without passive that squashes everything
figure(2) 
subplot(1,2,1) % for each mouse
y = [perc_modulated_loc; perc_modulated_replay]';
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
b(2).FaceColor = [0 0 0];

legend('tuning curve', 'halt replay')
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_modulated_loc; perc_modulated_replay]';
boxchart(y)
hold on

hold off
ylabel ('% modulated units')
xticklabels({'tuning curve', 'halt replay'})
legend('data','mean')

%% Putative sensorimotor Overlap of replay halt and location modulated
figure(3)
subplot(1,2,1)
bar(perc_overlap')
b(1).FaceColor = [0.5 0.5 0.5];
ylabel ('% modulated units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2)
boxchart(perc_overlap')
hold on
plot(mean(perc_overlap))
hold off
ylabel ('% modulated units')
xticklabels({'putative sensorimotor'})
