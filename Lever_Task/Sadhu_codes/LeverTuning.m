MyFilePath = 'D:\ImagingData\Cold Spots\Priyanka\Behavior\J4\J4_20200306_r0'

[Traces, TrialInfo,SingleUnits] = ParseBehaviorAndPhysiologyN(MyFilePath);

LeverBinSize = 0.25;
LeverBins = 0:LeverBinSize:5
odors = unique(TrialInfo.Odor);



K = load('D:\ImagingData\Cold Spots\Priyanka\Behavior\K1\processed\sing');
%% Tuning Curves - Lever 


 for o = 1:length(odors)   
     whichtrials = find(TrialInfo.Odor==odors(o));
for i = 1:length(K.SingleUnits)
    
    count = 0;
    for k = LeverBins
        count = count +1;
       
        count1 = 0;
        
    for j = 1:numel(whichtrials)
        
        
        if TrialInfo.TimeIndices(whichtrials(j),2) <= length(Traces.Lever{whichtrials(j)})
            k1 =  find(Traces.Lever{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):TrialInfo.TimeIndices(whichtrials(j),2),:) >= k & Traces.Lever{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):TrialInfo.TimeIndices(whichtrials(j),2),:) <(0.25+k));
        else
            k1 =  find(Traces.Lever{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):length(Traces.Lever{whichtrials(j)}),:) >= k & Traces.Lever{whichtrials(j)}(TrialInfo.TimeIndices(whichtrials(j),1):length(Traces.Lever{whichtrials(j)}),:) <(0.25+k));
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
    for k = LeverBins
        count = count +1;
        matr(count,1) = count;
        matr(count,2) = k;
    end
    


   for o = 1:numel(odors)
for i = 1:length(K.SingleUnits)
    for j = 1: numel(LeverBins)
        k = find(tunedspikespersec{o}{i}(j,:) ~= 0);
        if k~= 0
        k1 = tunedspikespersec{o}{i}(j,k);
        mottrials = numel(Traces.Lever) - nomotortrials{o}{i}(j);
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
   


save('C:\Users\sadhu\Documents\MATLAB\Analysis_Data\Priyanka\J4\Closed_Loop\TuningCurves_matrices\tuncurves_Lever','tunfinal','stdev');

%% Displaying Tuning Curves
for o = 1 : numel(odors)
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
    xlabel('Lever Volatge')
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
     xlabel('Lever Volatge')
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
     xlabel('Lever Volatge')
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
  
    
    
         