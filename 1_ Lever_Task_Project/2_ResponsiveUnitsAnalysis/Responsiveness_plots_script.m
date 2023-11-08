%Responsiveness_plots_script

%%
SessionName = {'O3/O3_20211005_r0_processed.mat',...
'O8/O8_20220702_r0_processed.mat',...
'O9/O9_20220630_r0_processed.mat',...
'S1/S1_20230314_r0_processed.mat',...
'S3/S3_20230321_r0_processed.mat',...
'S6/S6_20230718_r0_processed.mat',...
'S7/S7_20230707_r0_processed.mat',... %this session is a free lever one because the classic one is bugged - to be changed when fixed
'S11/S11_20230812_r0_processed.mat',...
'S12/S12_20230727_r0_processed.mat'};

for session = 1:length(SessionName)
    MouseName{session} = fileparts(SessionName{session});
    [responsive_units, resp_score] = GetResponsiveUnits(SessionName {session});
    
    Nb_resp_units(session) = sum(responsive_units);
    Nb_unit(session) = length(responsive_units);
    perc_responsive(session) = 100*Nb_resp_units(session)/Nb_unit(session);
    
    responsiveness.(['mouse',num2str(session)]).all_resp_score = resp_score;

    %getting response type
    odor_exc = 1*any(squeeze(resp_score)==1, [2 3]); %1 if excited by any event 
    odor_inh = 2*any(squeeze(resp_score)==2, [2 3]); %1 if inihibited by any event 
    exc_only = (odor_exc + odor_inh) == 1;
    inh_only = (odor_exc + odor_inh) == 2;
    mix  = (odor_exc + odor_inh) == 3;
    
    %getting odor specific responsiveness
    odor_resp = any(squeeze(resp_score)~=0, 3); %1 if showing a response for any event sorted by odor
    odor1_only = (odor_resp(:,1) -odor_resp(:,2) -odor_resp(:,3)) == 1;
    odor2_only = (odor_resp(:,2) -odor_resp(:,1) -odor_resp(:,3)) == 1;
    odor3_only = (odor_resp(:,3) -odor_resp(:,1) -odor_resp(:,2)) == 1;
    odor123_resp = (odor_resp(:,1) +odor_resp(:,2) +odor_resp(:,3)) ==3;
    odor12_resp = (odor_resp(:,1)+odor_resp(:,2)- odor_resp(:,3)) ==2;
    odor13_resp = (odor_resp(:,1)+odor_resp(:,3)- odor_resp(:,2)) ==2;
    odor23_resp = (odor_resp(:,2)+odor_resp(:,3)- odor_resp(:,1)) ==2;

    %getting event specific responsiveness
    event_resp = any(squeeze(resp_score)~=0, 2);
    trial_on = event_resp(:,1)==1;
    odor_on = event_resp(:,2)==1;
    trial_off = event_resp(:,3)==1;
    target = event_resp(:,5) + event_resp(:,end)~=0;

    %percentages
    %type of activation
    perc_exc(session) = 100*sum(exc_only)/Nb_unit(session);
    perc_inh(session) = 100*sum(inh_only)/Nb_unit(session);
    perc_mix(session) = 100*sum(mix)/Nb_unit(session);
    %odor specific
    perc_odor1(session) = 100*sum(odor1_only)/Nb_unit(session);
    perc_odor2(session) = 100*sum(odor2_only)/Nb_unit(session);
    perc_odor3(session) = 100*sum(odor3_only)/Nb_unit(session);
    perc_odor123(session) = 100*sum(odor123_resp)/Nb_unit(session);
    perc_odor12(session) = 100*sum(odor12_resp)/Nb_unit(session);
    perc_odor13(session) = 100*sum(odor13_resp)/Nb_unit(session);
    perc_odor23(session) = 100*sum(odor23_resp)/Nb_unit(session);
    %event specific
    perc_trialon(session) = 100*sum(trial_on)/Nb_unit(session);
    %perc_odoron - for now ignored as this response is often picked up by
    %trial on or target 
    perc_target(session) = 100*sum(target)/Nb_unit(session);
    perc_trialoff(session) = 100*sum(trial_off)/Nb_unit(session);
    
    %storing values for each mouse - currently not being used
    responsiveness.(['mouse',num2str(session)]).odor1 = odor1_only;
    responsiveness.(['mouse',num2str(session)]).odor2 = odor2_only;
    responsiveness.(['mouse',num2str(session)]).odor3 = odor3_only;
    responsiveness.(['mouse',num2str(session)]).odor123 = odor123_resp;
    responsiveness.(['mouse',num2str(session)]).odor12 = odor12_resp;
    responsiveness.(['mouse',num2str(session)]).odor13 = odor13_resp;
    responsiveness.(['mouse',num2str(session)]).odor23 = odor23_resp;
    responsiveness.(['mouse',num2str(session)]).excited = exc_only;
    responsiveness.(['mouse',num2str(session)]).inhibited = inh_only;
    responsiveness.(['mouse',num2str(session)]).mixed = mix;
    responsiveness.(['mouse',num2str(session)]).trialon = trial_on;
    responsiveness.(['mouse',num2str(session)]).trialoff = trial_off;
    responsiveness.(['mouse',num2str(session)]).target = target;

end 

%% Plots
figure(1) % Number of cells
boxchart(Nb_unit,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12,'XColor', 'none')
ylim([0 120])
yticks(20:20:120)
ylabel('Number of recorded units')
hold on
plot(mean(Nb_unit),'.','MarkerSize',12,'Color','k')
hold off
xticks([])

figure(2) % percentage responsive: excited + inhibited + mix
subplot(1,2,1) % for each mouse
y = [perc_exc; perc_inh; perc_mix]';
b = bar(y,'stacked');
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 1 1];
b(3).FaceColor = [0.5 0.5 0.5];

ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
legend('excited', 'inhibited', 'mixed')
ylabel ('% responsive units')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [perc_responsive; perc_exc ; perc_inh; perc_mix]';
boxchart(y,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% responsive units')
xticklabels({'responsive','excited','inhibited', 'mixed'})


figure(3)
y1 = [perc_odor1; perc_odor2; perc_odor3; perc_odor12+perc_odor13+perc_odor23; perc_odor123]';
y2 = [mean(perc_odor1); mean(perc_odor2); mean(perc_odor3); mean(perc_odor12+perc_odor13+perc_odor23); mean(perc_odor123)]';
boxchart(y1,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(y2,'.','MarkerSize',12,'Color','k')
xticklabels({'odor A only','odor B only','odor C only', '2 odors', '3 odors'})
ylabel ('% responsive units')

figure(4)
y1 = [perc_trialon; perc_target; perc_trialoff]';
y2 = [mean(perc_trialon); mean(perc_target); mean(perc_trialoff)]';
boxchart(y1,'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(y2,'.','MarkerSize',12,'Color','k')
xticklabels({'trial on','odor on','trial off'})
ylabel ('% responsive units')

figure(5) % mouse specific odor responses
for i = 1:length(MouseName)
    subplot(3, round(length(MouseName)/3),i)
    x = [perc_odor1(i) perc_odor2(i) perc_odor3(i) perc_odor12(i)+perc_odor13(i)+perc_odor23(i) perc_odor123(i)];
    labels = {'odor A only','odor B only','odor C only', '2 odors', '3 odors'} ;
    pie(x, labels)
end 

