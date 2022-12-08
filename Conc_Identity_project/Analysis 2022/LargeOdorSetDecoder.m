function LargeOdorSetDecoder(smoothPSTH)
% Discrimination across a larger (16) odor set: we trained classifiers to 
% perform binary odor discriminations of a target odorant from an 
% increasing number (1 to 15) of non-target odorant.

temp=smoothPSTH(:,:,:,1:5);
bins = 0:0.2:4;
Fs = 1000;

AONPre = [];
for k = 1:length(bins)-1
    start = Fs*(6 + bins(1)); %used to be 10 instead of 6
    stop = Fs*(6+bins(k+1)); %used to be 10 instead of 6

    AONPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));

end
CUMCUBEPre = AONPre;

%%

KernelType = 'Poly';
doPlot = 0;
nBoot = 3;
ANSWER_Pre = zeros(16,20,nBoot);
ANSWER_Mus = zeros(16,20,nBoot);

for odorNum = 2:16
    for boot = 1:nBoot
        odorSeq = randperm(16,odorNum);
        Results = SVM_OdorIdentifierLargeOdor(CUMCUBEPre(:,:,odorSeq,:),KernelType,doPlot);                
        ANSWER_Pre(odorNum-1,:,boot) = mean([Results.TruePos],2);  
    end
    odorNum
end
%% Plotting

figure(101);
subplot(2,1,1);hold on;
imagesc(squeeze(mean(ANSWER_Pre,3)),[0 1]);colormap('jet')
for j = 1:size(ANSWER_Pre,1)    
    figure(101);
    subplot(2,1,1);hold on;
    temp = squeeze(mean(ANSWER_Pre,3));
    idx = find(temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
        TP50(j,1) = idx(1);
    end
    axis([0 20 0 17])
     
end

end 

%% ---
%% MAIN FUNCTION
function [Results] = SVM_OdorIdentifierLargeOdor(CUMCUBEPre,Type,doPlot)

Results = struct('Params',{},'TruePos',{});
nRep = 4;
nStim = size(CUMCUBEPre,3);
timeSteps = 20;
for test = 1:nRep
    trainSet = setdiff(1:nRep,test);

    PERF_Pre = zeros(length(timeSteps),2);
    
    for t = 1:timeSteps
        
        x = squeeze(CUMCUBEPre(t,:,:,1:nRep));
        x = x(:,:);
         
        
        temp_x = reshape(x,[size(x,1) nStim nRep]);
        xTrain = temp_x(:,:,trainSet);
        xTest = temp_x(:,:,test);
        
        [~,ncategory,nrep] = size(xTrain);
        nstim = ncategory;

        Output = zeros(ncategory,nstim);

        for k = 1:ncategory
            Output(k,k:ncategory:end) = 1;
        end
        Output = repmat(Output,[1 nrep]);

        if doPlot
            figure;
        end
            tic
            XTrain = xTrain(:,:)';
            XTest = xTest';
            for k = 1:ncategory
                Y = Output(k,:)';

                SVMModel = fitcsvm(XTrain,Y,'Standardize',true,'KernelFunction',Type,...
                    'KernelScale','auto');
                [a,~]= predict(SVMModel,XTest(:,:));
                label_pre(k,:) = a;
                
            end

            perf_pre = mean(diag(label_pre));

            if (doPlot == 1)
                imagesc(label_pre,[0 1])
                colormap('cool')
                set(gcf,'WindowStyle','docked');
                title(['Performance is ' num2str(perf_pre) '%'])
                drawnow;
            end
            false_perf_pre = mean(~diag(label_pre));
            PERF_Pre(t,1) = perf_pre;
            PERF_Pre(t,2) = false_perf_pre;  
           
    end
    Results(test).TruePos = PERF_Pre(:,1);
    Results(test).MissClass = PERF_Pre(:,2);

end
end 