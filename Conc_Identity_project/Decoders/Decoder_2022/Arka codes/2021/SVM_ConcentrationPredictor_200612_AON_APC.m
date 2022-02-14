clear;clc;close all;
dataPath = 'D:\Dropbox\MATLAB\Olfaction\MarieData\210721_dataAON_APCexperiments\NewData\';

load([dataPath 'smoothPSTH200_AON.mat'])
temp=smoothPSTH(:,:,:,1:5);
bins = 0:0.2:4;
Fs = 1000;

AONPre = [];
for k = 1:length(bins)-1
    start = Fs*(10 + bins(1));
    stop = Fs*(10+bins(k+1));
%     start = Fs*(10+bins(k));
%     stop = start+200;
    AONPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
end
%% BASIC CONSISTENCY CHECKS and Preprocesing

% CUMCUBEPre = [];
% count = 1;
% temp = squeeze(sum(APCPre));
% for n = 1:size(APCPre,2)
%     flag=0;
%     for s = 1:20
%         idx = (temp(n,s,:)>0);
%         if (sum(idx)>=4)
%             tempidx = find(idx==1);
%         	CUMCUBEPre(:,count,s,1:4) = APCPre(:,n,s,tempidx(1:4));
%             flag = flag+1;
%         end
%     end
%     if (flag>=1)
%         count = count+1;
%     end
% end

CUMCUBEPre = [];
count = 1;
temp = squeeze(sum(AONPre));
for n = 1:size(AONPre,2)
    flag=0;
    for s = 1:20
        idx = (temp(n,s,:)>0);
        if (sum(idx)>=4)
            tempidx = find(idx==1);
        	CUMCUBEPre(:,count,s,1:4) = AONPre(:,n,s,tempidx(1:4));
            flag=flag+1;
        end
    end
    if (flag==20)
        count = count+1;
    end
end

%
figure;
temp = squeeze(sum(CUMCUBEPre));
for n = 1:25
    temp2 = (squeeze(temp(n,:,:))>0);
    subplot(5,5,n);imagesc(temp2)
end
colormap('gray')

%


% load([dataPath 'smoothPSTH200_APC.mat'])
% tempAPC=smoothPSTH(:,:,:,1:4);
% bins = 0:0.2:4;
% Fs = 1000;
% APCPre = [];
% for k = 1:length(bins)-1
%     start = Fs*(10 + bins(1));
%     stop = Fs*(10+bins(k+1));
% %     start = Fs*(10+bins(k));
% %     stop = start+200;
%     APCPre(k,:,:,:) = squeeze(nanmean(tempAPC(:,:,start:stop,:),3));
% end
% CUMCUBEPre = APCPre;

%% BASIC PLOTTING
figure;
for k = 1:20
    subplot(5,4,k)
    temp = squeeze(CUMCUBEPre(:,200+k,12,:));
    plot(temp,'LineWidth',2);
%     line([40 40],[0 max(temp(:))]); hold on;
%     line([60 60],[0 max(temp(:))]);
%     axis([0 65 0 1+max(temp(:))])
end
%
% figure;
% for k = 1:20
%     subplot(5,4,k)
%     temp = squeeze(APCPre(:,200+k,11,:));
%     plot(temp,'LineWidth',2);
% %     line([5 5],[0 max(temp(:))]); hold on;
% %     line([25 25],[0 max(temp(:))]);
%     axis([0 20 0 1+max(temp(:))])
% end
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
                Xtest_mus = xtest_mus(CellId,:)';
                for k = 1:ncategory
                    Y = Output(k,:)';
                    template = templateSVM('Standardize',true);
                    
                    SVMModel = fitcecoc(X,Y,'Learners',template);
 

                    [a,~]= predict(SVMModel,Xtest(:,:));
                    label(k,:) = a;

                    [a,~]= predict(SVMModel,Xtest_mus(:,:));
                    label_mus(k,:) = a;

                end
                PERF(boot,1) = corr2(label,Output(:,1:20));
                PERF(boot,2) = corr2(label_mus,Output(:,1:20));
            end
%             GENERAL_PERF_Pre(t,CellIter,test) = mean(PERF(:,1));
%             GENERAL_PERF_Mus(t,CellIter,test) = mean(PERF(:,2));
            GENERAL_PERF_Pre(t,CellIter,test,:) = squeeze(PERF(:,1));
            GENERAL_PERF_Mus(t,CellIter,test,:) = squeeze(PERF(:,2));
        end
    end
    t
end
%% PLOTTING FOR DECODERS
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*nboot]);
temp2 = reshape(GENERAL_PERF_Mus,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*nboot]);

cellNumIdx = 42;
shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp1(:,cellNumIdx,:)),2),...
    nanstd(squeeze(temp1(:,cellNumIdx,:)),0,2)./sqrt(4*nboot),'LineProps','--k')
shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp2(:,40,:)),2),...
    nanstd(squeeze(temp2(:,cellNumIdx,:)),0,2)./sqrt(4*nboot),'LineProps','--r')

%
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[20 size(CellNum,2) 4*nboot]);
temp2 = reshape(GENERAL_PERF_Mus,[20 size(CellNum,2) 4*nboot]);

TimeBinIdx = 8;
shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp1(TimeBinIdx,:,:)),2),...
    nanstd(squeeze(temp1(TimeBinIdx,:,:)),0,2)./sqrt(12),'LineProps','-k')
shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp2(TimeBinIdx,:,:)),2),...
    nanstd(squeeze(temp2(TimeBinIdx,:,:)),0,2)./sqrt(12),'LineProps','--r')

figure;imagesc(squeeze(nanmean(temp1,3))',[0 1]);colormap('jet')