%WaveformAmplitudeViewPlots
%% paths
SessionName = fullfile('S12','S12_20230727_r0_processed.mat');

if strcmp(computer, 'MACI64')
    ProcessedBehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/';
else
    ProcessedBehaviorPath = '/mnt/data/Processed/Behavior/';
end

MySession = fullfile(ProcessedBehaviorPath,SessionName);

%% Get the data loaded
[TracesOut, ColNames, handles.TrialInfo, SingleUnits, TTLs, ...
    ReplayTTLs, SampleRate, TimestampAdjuster, PassiveTracesOut, StartStopIdx, OpenLoop] = ...
    LoadProcessedDataSession(MySession); % LoadProcessedSession; % loads relevant variables

handles.SingleUnits = SingleUnits;
handles.TimestampAdjuster = TimestampAdjuster;
Nb_unit = size(SingleUnits,2);



%% Plot spike amplitudes
units_per_fig = 5;
for whichUnit = 1:Nb_unit % for every cell
            if mod(whichUnit,units_per_fig)
                FRplot = 2*mod(whichUnit,units_per_fig);
            else
                FRplot = 2*units_per_fig;
            end
            Rasterplot = FRplot - 1;

subplot(units_per_fig,2,Rasterplot);
thisunitamps = handles.SingleUnits(1).spikescaling(find(handles.SingleUnits(1).clusterscalingorder == handles.SingleUnits(whichUnit).id));
%axes(handles.amplitudeaxes);
plot(handles.SingleUnits(whichUnit).spikes,thisunitamps,'.k');
ylim([0 50])
set(gca, 'TickDir', 'out')
hold on
session_end = handles.TrialInfo.SessionTimestamps(end,2) + handles.TimestampAdjuster;
line([session_end session_end],get(gca,'YLim'),'Color','k');
hold off

% Get tetrode and cluster id 
TetrodeOrder(whichUnit,:) = [SingleUnits(whichUnit).tetrode SingleUnits(whichUnit).id]; 
handles.tetrode.String = num2str(TetrodeOrder(whichUnit,1));
handles.Cluster_ID.String = num2str(TetrodeOrder(whichUnit,2));
title(['Unit# ',num2str(whichUnit), '; Clust# ', num2str(SingleUnits(whichUnit).id),...
                            '; tet# ',num2str(SingleUnits(whichUnit).tetrode)]);

%Get the amplitude distributions for behavior and passive parts of the
%session
subplot(units_per_fig,2,FRplot);
i = find(handles.SingleUnits(whichUnit).spikes > session_end, 1); % index of 1st spike after behavior ends
histogram(thisunitamps(1:i-1),'BinWidth', 0.5,'BinLimits',[0,50],'Normalization', 'probability')
hold on
histogram (thisunitamps(i:end),'BinWidth', 0.5,'BinLimits',[0,50],'Normalization', 'probability')

if mod(whichUnit,units_per_fig) == 0
    figure;
end

end