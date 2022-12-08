function OdorGeneralizer(smoothPSTH)
% Generalization to a novel concentration: after odor identity was learned 
% from a test set of concentration, we trained the classifier to identify 
% the odorant (one of a total five possibilities) by learning on three out 
% the four concentrations, and further tested for generalization to a novel
% held-out concentration.

temp=smoothPSTH(:,:,:,:,3:7);
temp = reshape(temp,[size(smoothPSTH,1),20,10000,5]);
bins = 0:0.2:4;
Fs = 1000;

CUMCUBEPre = [];
for k = 1:length(bins)-1
    start = Fs*(6 + bins(1)); %used to be 10 instead of 6
    stop = Fs*(6+bins(k+1)); %used to be 10 instead of 6

    CUMCUBEPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
end

%% DECODER

ncells = size(CUMCUBEPre,2);
CellNum = 10:5:ncells;

nboot = 10;
GENERAL_PERF_Pre = zeros(20,length(CellNum),4,nboot,5);

for t = 1:size(CUMCUBEPre,1)
    t
    RESPCUBE_Pre = squeeze(CUMCUBEPre(t,:,:,1:4));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 4]); 
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 4]);
       
    for test = 1:4
        trainset = setdiff(1:4,test);
        xtrain = squeeze(x(:,:,trainset,:));

        xtest = squeeze(x(:,:,test,:));
        
        [nroi,ncategory,ndil,nrep] = size(xtrain);
        nstim = ncategory*ndil;
        nTrainRep = nrep;

        Output = zeros(ncategory,nstim);
        Output_cross = zeros(ncategory,nstim,nrep);

        for k = 1:ncategory
            Output(k,k:ncategory:end) = 1;
        end
        Output = repmat(Output,[1 nrep]);
                
        for CellIter = 1:length(CellNum)
            PERF = zeros(nboot,2,5);            
            for boot = 1:nboot
                CellId = randperm(ncells,CellNum(CellIter));
                tic
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';
                
                for k = 1:ncategory
                    Y = Output(k,:)';
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto');

                    [a,~]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[ncategory,nrep]);
                    label(k,:) = mean(temp,2);

                end
                PERF(boot,1,1:5) = diag(label);
                GENERAL_PERF_Pre(t,CellIter,test,boot,1:5) = squeeze(PERF(boot,1,:));

            end

        end

    end
end
%% PLOTTING FOR DECODERS

AONperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));

xx = squeeze(nanmean(GENERAL_PERF_Pre,4));

figure;
figCount = 1;
for j = 1:5
    subplot(3,2,figCount);
    imagesc(squeeze(mean(xx(:,:,:,j),3))',[0 1])
    axis('square');
    figCount = figCount+1;
    xlim([2.5 20.5])
    set(gca,'TickDir','out')
end
colormap('jet')

figure;
imagesc(squeeze(mean(xx(:,:,:,:),[3 4]))',[0 1])
axis('square');
xlim([2.5 20.5])
set(gca,'TickDir','out')
colormap('jet')
end 