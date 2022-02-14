function [Results] = SVM_OdorIdentifierLargeOdor_190908_AON_APC(CUMCUBEPre,Type,doPlot)

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
        
       