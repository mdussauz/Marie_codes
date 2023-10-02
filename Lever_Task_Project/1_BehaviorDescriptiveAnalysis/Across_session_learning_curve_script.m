%% File path
MouseName = 'S3';
BehaviorPath = '/Users/mariedussauze/Desktop/Analysis/data/Behavior';
MyFilePath = fullfile(BehaviorPath,MouseName);

%% Getting variables
[NumTrials, SuccessRate, MaxTargetHold, TotalTargetStay,  ExpertReached,...
    AllMotorLocInTrial, AllCenteredLeverInTrial]...
    = GetLearningCurveInfo(MyFilePath);


%% Plotting

fig2 = figure(2); %histograms of motor location in trials for all sessions
BinWidth = 8;
bins = -100:BinWidth:100; %Bin width of 8 allows to have entire tz for motor as one of the bin
for j = 1: ExpertReached
    subplot(ceil(ExpertReached/5),5,j, 'Parent', fig2); hold on
    if any(~isnan(AllMotorLocInTrial{1,j}))
        [N,edges] = histcounts(AllMotorLocInTrial{1,j}, bins);
        histogram('BinEdges',edges,'BinCounts',N/sum(N))
    end 
end
hold off
%NB: Motor histogram cannot be used for further analysis as some sessions
%had failings of the rotary encoder

fig3 = figure(3); % histograms of lever location in trials for all sessions
BinWidth = 0.6;
bins = -4.5:BinWidth:4.5;% Bin width of .6 allows to have entire tz for lever as one of the bins
%TZCenterInd = find(BinCenters==0);
BinCenters = (bins(1:end-1)+bins(2:end))/2;
fwhmx = NaN(ExpertReached,1);
for j = 1: ExpertReached
    subplot(ceil(ExpertReached/5),5,j, 'Parent', fig3); hold on
    if any(~isnan(AllCenteredLeverInTrial{1,j}))
        [N,edges] = histcounts(AllCenteredLeverInTrial{1,j}, bins);
        %h = histogram('BinEdges',edges,'BinCounts',(N/sum(N))*100);
        occupancy = N/sum(N);
        plot(BinCenters,(occupancy*100))

        % -- Calculate full width at half occupancy in target zone (i.e. -0.3
        % to 0.3) 
        % 1) extrapolate to get more points - otherwise might get the tz
        % occupancy value instead of half of that
        extrapt=10;
        xq = linspace(BinCenters(1),BinCenters(end),length(BinCenters)*extrapt);
        yq = interp1(BinCenters,occupancy,xq);
        plot(xq,yq*100)
        % 2) half occupancy at tz
        TzOccupancy = yq(75)/2;
        % 3) Find left index of half the occupancy at tz
        index1 = find(yq >= TzOccupancy, 1, 'first');
        % 4) Find right index of half the occupancy at tz
        index2 = find(yq >= TzOccupancy, 1, 'last');
        % 5) Full width half max divided by total spread of lever position
        fwhmx(j) = (xq(index2) - xq(index1))/(BinCenters(end)-BinCenters(1));
        

    end
end
hold off

fig1 = figure(1); % Learning curve plots
subplot(3,2,1); hold on; % number of trials per session
plot(NumTrials)
ylabel('nb of trials')

subplot(3,2,2); hold on; % success rate per session
plot(SuccessRate)
ylabel('success rate (%)')

subplot(3,2,3); hold on %max time spent consecutively in tz
plot(MaxTargetHold)
ylabel('Mean of max target hold (ms)')

subplot(3,2,4); hold on %total time spent consecutively in tz
plot(TotalTargetStay)
ylabel('Mean of tot target stay (ms)')

subplot(3,2,5)
plot(1 - fwhmx)
ylabel('full width at half tz occupancy')



% figure settings
for i = 1:5
    subplot(3,2,i)
    xlabel('session number')
    set(gca,'TickDir','out');
end 
hold off




