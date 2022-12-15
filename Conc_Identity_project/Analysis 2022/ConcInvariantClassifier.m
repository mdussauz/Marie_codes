function [GENERAL_PERF_Pre] = ConcInvariantClassifier(smoothPSTH)
% Concentration invariant odor recognition: we trained the classifier to 
% assign odor identity irrespective of concentration. 
%%
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
Results = struct('Params',{},'OutputFit',{},'OutputCross',{}, ...
    'WeightsFit',{}, 'WeightsCross', {});

ncells = size(CUMCUBEPre,2);
CellNum = 10:5:ncells;


GENERAL_PERF_Pre = zeros(20,length(CellNum),4,5);
GENERAL_PERF_Mus = zeros(20,length(CellNum),4,5);

for t = 1:20
    t
    RESPCUBE_Pre = squeeze(CUMCUBEPre(t,:,:,1:4));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 4]);
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 4]);
    
    RESPCUBE_Mus = squeeze(CUMCUBEPre(t,:,:,1:4));
    x_mus = reshape(RESPCUBE_Mus,[ncells 5 4 4]);
    x_mus = zscore(x_mus(:,:),0,2);
    x_mus = reshape(x_mus,[ncells 5 4 4]);
    
    for test = 1:4
        trainset = setdiff(1:4,test);
        xtrain = squeeze(x(:,:,:,trainset));
        
        xtest = squeeze(x(:,:,:,test));
        xtest_mus = squeeze(x_mus(:,:,:,test));
        
        [nroi,ncategory,ndil,nrep] = size(xtrain);
        nstim = ncategory*ndil;
        nTrainRep = nrep;
        nboot = 10;
        
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
                Xtest_mus = xtest_mus(CellId,:)';
                for k = 1:ncategory
                    Y = Output(k,:)';
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto');
                    
                    [a,b]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[ncategory,ndil]);
                    label(k,:) = mean(temp,2);
                    
                    [a,b]= predict(SVMModel,Xtest_mus(:,:));
                    temp = reshape(a,[ncategory,ndil]);
                    label_mus(k,:) = mean(temp,2);
                    
                end
                PERF(boot,1,1:5) = diag(label);
                
                PERF(boot,2,1:5) = diag(label_mus);
            end
            GENERAL_PERF_Pre(t,CellIter,test,1:5) = mean(PERF(:,1,:));
            GENERAL_PERF_Mus(t,CellIter,test,1:5) = mean(PERF(:,2,:));
            %             CellNum(CellIter)
        end
        
    end
end

%% ----- plotting
figure;
figCount = 1;
for j = 1:5
    subplot(5,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Pre(:,:,:,j),3))',[0 1]) %mean across test
    %     axis('square');
    figCount = figCount+1;
%     subplot(5,2,figCount);
%     imagesc(squeeze(mean(GENERAL_PERF_Mus(:,:,:,j),3))',[0 1])
%     %     axis('square');
%     figCount = figCount+1;
end
colormap('jet')

figure; 
imagesc(squeeze(mean(GENERAL_PERF_Pre(:,:,:,1:4),[3 4]))',[0 1]) %mean across test and odor
colormap('jet')

end