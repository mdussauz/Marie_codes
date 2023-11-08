% clear;clc;close alL;
% dataPath = 'D:\Dropbox\MATLAB\Olfaction\MarieData\210721_dataAON_APCexperiments\';
% load([dataPath 'smoothPSTH200_AON.mat'])
% tempAON=smoothPSTH(:,:,:,1:5);
% bins = -1:0.2:5;
% Fs = 1000;
% 
% AONPre = [];
% for k = 1:length(bins)-1
% %     start = Fs*(10 + bins(1));
% %     stop = Fs*(10+bins(k+1));
%     start = Fs*(10+bins(k));
%     stop = start+200;
%     AONPre(k,:,:,:) = squeeze(nanmean(tempAON(:,:,start:stop,:),3));
% end
% CUMCUBEPre = AONPre;
% 
% load([dataPath 'smoothPSTH200_APC.mat'])
% tempAPC=smoothPSTH(:,:,:,1:5);
% 
% APCPre = [];
% for k = 1:length(bins)-1
% %     start = Fs*(10 + bins(1));
% %     stop = Fs*(10+bins(k+1));
%     start = Fs*(10+bins(k));
%     stop = start+200;
%     APCPre(k,:,:,:) = squeeze(nanmean(tempAPC(:,:,start:stop,:),3));
% end
% CUMCUBEPre = APCPre;
%%
clear;clc;close all;
dataPath = 'D:\Dropbox\MATLAB\Olfaction\MarieData\210721_dataAON_APCexperiments\NewData\';

load([dataPath 'smoothPSTH200_APC.mat'])
temp=smoothPSTH(:,:,:,1:5);
bins = 0:0.2:4;
Fs = 1000;

APCPre = [];
for k = 1:length(bins)-1
    start = Fs*(10 + bins(1));
    stop = Fs*(10+bins(k+1));
%     start = Fs*(10+bins(k));
%     stop = start+200;
    APCPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
end
%% BASIC PREPROCESSING STEPS

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
%             flag=flag+1;
%         end
%     end
%     if (flag==20)
%         count = count+1;
%     end
% end

CUMCUBEPre = APCPre;
%% DECODER

ncells = size(CUMCUBEPre,2);
CellNum = 10:5:ncells;

nboot = 10;
GENERAL_PERF_Pre = zeros(20,length(CellNum),4,nboot,5);
GENERAL_PERF_Mus = zeros(20,length(CellNum),4,nboot,5);

for t = 1:size(CUMCUBEPre,1)
    t
    RESPCUBE_Pre = squeeze(CUMCUBEPre(t,:,:,1:4));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 4]); 
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 4]);
    
    RESPCUBE_Mus = squeeze(CUMCUBEPre(t,:,:,1:4));
    x_mus = reshape(RESPCUBE_Pre,[ncells 5 4 4]); 
    x_mus = zscore(x_mus(:,:),0,2);
    x_mus = reshape(x_mus,[ncells 5 4 4]);
    
    for test = 1:4
        trainset = setdiff(1:4,test);
        xtrain = squeeze(x(:,:,trainset,:));

        xtest = squeeze(x(:,:,test,:));
        xtest_mus = squeeze(x_mus(:,:,test,:));
        
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
                Xtest_mus = xtest_mus(CellId,:)';
                for k = 1:ncategory
                    Y = Output(k,:)';
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto');

                    [a,~]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[ncategory,nrep]);
                    label(k,:) = mean(temp,2);

                    [a,b]= predict(SVMModel,Xtest_mus(:,:));
                    temp = reshape(a,[ncategory,nrep]);
                    label_mus(k,:) = mean(temp,2);

                end
                PERF(boot,1,1:5) = diag(label);
                PERF(boot,2,1:5) = diag(label_mus);
                GENERAL_PERF_Pre(t,CellIter,test,boot,1:5) = squeeze(PERF(boot,1,:));
                GENERAL_PERF_Mus(t,CellIter,test,boot,1:5) = squeeze(PERF(boot,2,:));
            end

        end

    end
end
%% PLOTTING FOR DECODERS

AONperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));
AONperf_mus = squeeze(nanmean(GENERAL_PERF_Mus,3));

xx = squeeze(nanmean(GENERAL_PERF_Pre,4));
yy = squeeze(nanmean(GENERAL_PERF_Mus,4));
figure;
figCount = 1;
for j = 1:5
    subplot(5,2,figCount);
    imagesc(squeeze(mean(xx(:,:,:,j),3))',[0 1])
    axis('square');
    figCount = figCount+1;
    subplot(5,2,figCount);
    imagesc(squeeze(mean(yy(:,:,:,j),3))',[0 1])
    axis('square');
    figCount = figCount+1;
    xlim([2.5 20.5])
    set(gca,'TickDir','out')
end
colormap('jet')


%
% figure;hold on;
% CellIdx = 10;
% temp1 = squeeze(AONperf_pre(:,CellIdx,:,4));
% shadedErrorBar(1:size(temp1,1),mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')

%
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*5*nboot]);
temp2 = reshape(GENERAL_PERF_Mus,[size(GENERAL_PERF_Pre,1) size(CellNum,2) 4*5*nboot]);

cellNumIdx = 42;
shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp1(:,cellNumIdx,:)),2),...
    nanstd(squeeze(temp1(:,cellNumIdx,:)),0,2)./sqrt(4*5*nboot),'LineProps','--k')
shadedErrorBar(1:size(GENERAL_PERF_Pre,1),nanmean(squeeze(temp2(:,cellNumIdx,:)),2),...
    nanstd(squeeze(temp2(:,cellNumIdx,:)),0,2)./sqrt(4*5*nboot),'LineProps','--r')
%
figure;hold on;
TimeBinIdx = 8;
shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp1(TimeBinIdx,:,:)),2),...
    nanstd(squeeze(temp1(TimeBinIdx,:,:)),0,2)./sqrt(4*5*nboot),'LineProps','-k')
shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp2(TimeBinIdx,:,:)),2),...
    nanstd(squeeze(temp2(TimeBinIdx,:,:)),0,2)./sqrt(4*5*nboot),'LineProps','--r')

figure;imagesc(squeeze(nanmean(temp1,3))',[0 1]);colormap('jet')
