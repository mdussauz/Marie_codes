function [GENERAL_PERF_Pre] = ConcentrationPredictor(smoothPSTH)
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

ncells = size(CUMCUBEPre,2);
CellNum = 10:5:ncells;

nboot = 10;

GENERAL_PERF_Pre = zeros(size(CUMCUBEPre,1),length(CellNum),4,nboot);
GENERAL_PERF_Mus = zeros(size(CUMCUBEPre,1),length(CellNum),4,nboot);

for t = 1:size(CUMCUBEPre,1)
    t
    RESPCUBE_Pre = squeeze(CUMCUBEPre(t,:,:,1:4));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 4]); 
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 4]);
       
    for test = 1:4
        trainset = setdiff(1:4,test);
        xtrain = squeeze(x(:,:,:,trainset));

        xtest = squeeze(x(:,:,:,test));

        
        [nroi,ncategory,ndil,nrep] = size(xtrain);
        nstim = ncategory*ndil;
        nTrainRep = nrep;
 
        Output = zeros(ncategory,nstim);
        Output_cross = zeros(ncategory,nstim,nrep);
 
        for k = 1:ncategory
            Output(k,k:ncategory:end) = [1 2 3 4];
        end
        Output = repmat(Output,[1 nTrainRep]);

        for CellIter = 1:length(CellNum)
            PERF = zeros(nboot,2);            
            for boot = 1:nboot
                CellId = randperm(ncells,CellNum(CellIter));
                tic
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';

                for k = 1:ncategory
                    Y = Output(k,:)';
                    template = templateSVM('Standardize',true);
                    
                    SVMModel = fitcecoc(X,Y,'Learners',template);
 

                    [a,~]= predict(SVMModel,Xtest(:,:));
                    label(k,:) = a;

                end
                PERF(boot,1) = corr2(label,Output(:,1:20));

            end

            GENERAL_PERF_Pre(t,CellIter,test,:) = squeeze(PERF(:,1));

        end
    end
    t
end
%% PLOTTING FOR DECODERS
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*nboot]);
%temp2 = reshape(GENERAL_PERF_Mus,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*nboot]);

cellNumIdx = 42;
shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp1(:,cellNumIdx,:)),2),...
    nanstd(squeeze(temp1(:,cellNumIdx,:)),0,2)./sqrt(4*nboot),'LineProps','--k')
% shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp2(:,40,:)),2),...
%     nanstd(squeeze(temp2(:,cellNumIdx,:)),0,2)./sqrt(4*nboot),'LineProps','--r')

%
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[20 size(CellNum,2) 4*nboot]);
%temp2 = reshape(GENERAL_PERF_Mus,[20 size(CellNum,2) 4*nboot]);

TimeBinIdx = 8;
shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp1(TimeBinIdx,:,:)),2),...
    nanstd(squeeze(temp1(TimeBinIdx,:,:)),0,2)./sqrt(12),'LineProps','-k')
% shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp2(TimeBinIdx,:,:)),2),...
%     nanstd(squeeze(temp2(TimeBinIdx,:,:)),0,2)./sqrt(12),'LineProps','--r')

figure;imagesc(squeeze(nanmean(temp1,3))',[0 1]);colormap('jet')

end 