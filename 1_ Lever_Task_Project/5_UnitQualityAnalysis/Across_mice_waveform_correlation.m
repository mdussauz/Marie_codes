%Across_mice_waveform_correlation_script

%% Classic open loop sessions
SessionName = {...
%'O3/2021-10-05_14-24-31',...
%'O8/2022-07-02_14-00-28',...
%'O9/2022-06-30_15-14-32',...
'S1/2023-03-14_15-01-07',...
'S3/2023-03-29_16-23-54',... %this is a mid session OL
'S6/2023-07-27_14-16-58',... 
'S7/2023-07-07_13-58-44',... 
'S11/2023-08-05_15-16-13',... %this is a free lever OL
'S12/2023-07-27_17-18-37'};

for session = 1:length(SessionName)
    mouse = session;
    MouseName{session} = fileparts(SessionName{session});

    [WFcorrelation, concWFcorrelation, CorrelationCriteria, concCorrelationCriteria] = GetWaveformCorrelation (SessionName{session});

    nb_units_total(mouse) = length(concWFcorrelation);
    %based on concatenated waveforn across the channels of the unit's main
    %tetrode:
    bad_units_conc = find(concWFcorrelation > concCorrelationCriteria);
    nb_units_conc(mouse) = sum(concWFcorrelation > concCorrelationCriteria);

    %based on single waveform on the unit's main channel:
    bad_units_single = find(WFcorrelation > CorrelationCriteria);
    nb_units_single(mouse) = sum(WFcorrelation > CorrelationCriteria);
end

%% Plot

%% for concatenated waveforms
figure(1) 
subplot(1,2,1) % for each mouse
y = [100*nb_units_conc./nb_units_total];
b = bar(y);
b(1).FaceColor = [0.5 0.5 0.5];
ylim([0 100])
yticks(20:20:100)
set(gca,'box','off','color','none','TickDir','out','linewidth',2,...
    'fontname','calibri','fontsize',12)
ylabel ('% units with stable waveform')
xticklabels(MouseName)
xlabel('Mouse Name')

subplot(1,2,2) % across mice
y = [100*nb_units_conc./nb_units_total]';
boxchart(y, 'MarkerColor',"#A2142F", 'BoxWidth',0.4, 'BoxFaceColor','k');
set(gca,'box','off','color','none','TickDir','out','XTickLabelRotation' ,45,'linewidth',2,...
    'fontname','calibri','fontsize',12)
ylim([0 100])
yticks(20:20:100)
hold on
plot(mean(y),'.','MarkerSize',12,'Color','k')
hold off
ylabel ('% units with stable waveform')

