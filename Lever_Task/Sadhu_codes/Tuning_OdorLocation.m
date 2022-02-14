MyFilePath = '/Users/mdussauz/Desktop/Analysis/M04/2021-03-30_16-46-03/M04_20210330_r0';
MyFileName = 'M04_20210330_r0';
myephysdir = '/Users/mdussauz/Desktop/Analysis/M04/2021-03-30_16-46-03';
[Traces,TrialInfo] = ParseBehaviorAndPhysiology(MyFilePath);
%ParseBehaviorAndPhysiology(MyBehaviorFile);

%% core data extraction (and settings)
[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

sessionstart = MyData(1,1);
sessionstop = MyData(end,1);

%% Parse into trials
[Trials] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[Traces, TrialInfo, TargetZones] = ParseBehaviorTrials(MyData, MySettings, DataTags, Trials, sessionstart, sessionstop);

%% Get info from the OEPS files if available
%[myephysdir] = WhereSpikeFile(MyFileName);

% get all TTLs for the open ephys session

[~,TTLs] = GetOepsAuxChannels(myephysdir, Trials.TimeStamps, 'ADC', 0);


%% Get spikes - label spikes by trials

SingleUnits = GetSingleUnits(myephysdir);
[SingleUnits] = Spikes2Trials(TTLs, SingleUnits);
%[SingleUnits, EphysTuningTrials] = Spikes2Trials_Tuning(myephysdir, TS, TrialInfo, MyTuningTrials);

% keep only good units
GoodUnits = [];
for i = 1:size(SingleUnits,2)
    if SingleUnits(i).quality == 1
        GoodUnits = [GoodUnits i];
    end
end


%%
save('/Users/mdussauz/Desktop/Analysis/M04/processed/sing','SingleUnits');

MotorBinSize = 10;
MotorBins = -100:MotorBinSize:100
odors = unique(TrialInfo.Odor);

K = load('/Users/mdussauz/Desktop/Analysis/M04/processed/sing');

%% Tuning Curves - Odor Location
for o = 1:numel(odors)
     whichtrials = find(TrialInfo.Odor==odors(o));
for i = 1:length(K.SingleUnits)
    
    count = 0;
    for k = MotorBins
        count = count +1;
       
        count1 = 0;
    for j = 1:numel(whichtrials)
        
        
        if TrialInfo.TimeIndices(whichtrials(j),2) <= length(Traces.Motor{whichtrials(j)})
            k1 =  find(Traces.Motor{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):TrialInfo.TimeIndices(whichtrials(j),2),:) >= k & Traces.Motor{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):TrialInfo.TimeIndices(whichtrials(j),2),:) <(MotorBinSize+k));
        else
            k1 =  find(Traces.Motor{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):length(Traces.Motor{whichtrials(j)}),:) >= k & Traces.Motor{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):length(Traces.Motor{whichtrials(j)}),:) <(MotorBinSize+k));
        end
            k1 = k1*0.002;
            s(:,1) = K.SingleUnits(i).trialalignedspikes;
            s(:,2) = K.SingleUnits(i).trialtags;
            tsp = find(s(:,2) == whichtrials(j));
            trialspi = s(tsp,1);
            spi = 0;
            if ~isempty(k1)
            for j1 = 1:length(k1)
                for j11 = 1:length(trialspi)
                    if abs(k1(j1) - trialspi(j11)) <= 0.002
                        spi = spi+1;
                    end
                end
            end
            
            
            spikespersec = spi/(0.002*length(k1));
            else
                spikespersec = 0;
                count1 = count1 + 1;
                
                
            end
            tunedspikespersec{o}{i}(count,j) = spikespersec;
            tunedspikes{o}{i}(count,j) = spi;
            
    end
    nomotortrials{o}{i}(count) = count1;
    clear k1;
        clear s;
        clear tsp;
    end
end
end

count = 0;
    for k = -100:10:100
        count = count +1;
        matr(count,1) = count;
        matr(count,2) = k;
    end
    
    


 for o = 1:numel(odors)  
for i = 1:length(K.SingleUnits)
    for j = 1 : numel(MotorBins)
        k = find(tunedspikespersec{o}{i}(j,:) ~= 0);
        if k~= 0
        k1 = tunedspikespersec{o}{i}(j,k);
        mottrials = numel(TrialInfo.TrialID) - nomotortrials{o}{i}(j);
        k11 = zeros(1,mottrials);
        k11(1:length(k1)) = k1(:);
        
        tunfinal{o}{i}(j) = mean(k11);
        stdev{o}{i}(j) = std(k11)/sqrt(length(k11));
        else
            tunfinal{o}{i}(j) = 0;
            stdev{o}{i}(j) = 0;
        end
        
        
    end
    
    clear k
    
    clear k11
   
end
 end


save('/Users/mdussauz/Desktop/Analysis/M04/TuningCurves_matrices/tuncurves_OdorMotor','tunfinal','stdev');

%% Displaying Tuning Curves
for o = 1:numel(odors)
for i = 1:length(K.SingleUnits)
    if i<= 25
    subplot(5,5,i)
    f1 = figure(1)
    set(f1,'Position',[1 1 1680 1050]);
     set(f1, 'visible', 'off');
     shadedErrorBar(matr(:,2),tunfinal{o}{i},stdev{o}{i},'lineProps','r')
     hold on
    plot(matr(:,2),tunfinal{o}{i},'Color','black','LineWidth',1.2')
    set(gcf,'color','w');
    xlabel('Odor Motor Position')
    ylabel('Spikes/sec')
    k = strcat('Unit :',{' '},num2str(i));
    title(k)
    end
    
    if i>25 && i<=50
        f2 = figure(2)
        set(f2,'Position',[1 1 1680 1050]);
         set(f2, 'visible', 'off');
        subplot(5,5,i-25)
     shadedErrorBar(matr(:,2),tunfinal{o}{i},stdev{o}{i},'lineProps','r')
     hold on
    plot(matr(:,2),tunfinal{o}{i},'Color','black','LineWidth',1.2')
    set(gcf,'color','w');
    xlabel('Odor Motor Position')
    ylabel('Spikes/sec')
    k = strcat('Unit :',{' '},num2str(i));
    title(k)
    end
    
  if i > 50
         f3 = figure(3)
        set(f3,'Position',[1 1 1680 1050]);
         set(f3, 'visible', 'off');
        subplot(5,5,i-50)
    shadedErrorBar(matr(:,2),tunfinal{o}{i},stdev{o}{i},'lineProps','r')
     hold on
    plot(matr(:,2),tunfinal{o}{i},'Color','black','LineWidth',1.2')
    set(gcf,'color','w');
    xlabel('Odor Motor Position')
    ylabel('Spikes/sec')
    k = strcat('Unit :',{' '},num2str(i));
    title(k)
    
        
        
        
    
  
    
  end
end
       saveas(f1,strcat('Tun_Odor_',num2str(o),'_1'),'png');
    saveas(f2,strcat('Tun_Odor_',num2str(o),'_2'),'png');
%     saveas(f3,strcat('Tun_Odor_',num2str(o),'_3'),'png');
    close(f1)
    close(f2)
%     close(f3)
end
   
        
  