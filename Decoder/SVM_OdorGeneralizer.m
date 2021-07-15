clear;clc;close aLL;
% dataPath = 'D:\Dropbox\MATLAB\ChaeBanerjeeAnalysis\NatureSubmission\Data\Concentration\';
% load([dataPath 'MitralCellsPreMuscimol.mat'])
% % permuteIdx = randperm(20);
% 
% CUMCUBEPre=CumulativeCube(:,:,:,1:4);
% 
% load([dataPath 'MitralCellsMuscimol.mat'])
% CUMCUBEMus=CumulativeCube(:,:,:,1:4);
% 
% Results = struct('Params',{},'OutputFit',{},'OutputCross',{}, ...
%     'WeightsFit',{}, 'WeightsCross', {});
%---------------------------------------------------------------------------

dataPath = 'D:\Dropbox\MATLAB\ChaeBanerjeeAnalysis\NatureSubmission\Data\Concentration\';
load([dataPath 'MitralCellsPreMuscimol.mat'])
CUMCUBEPre=CumulativeCube(:,:,:,1:4);
dataPath = 'D:\Dropbox\MATLAB\ChaeBanerjeeAnalysis\NatureSubmission\Data\Concentration\';
load([dataPath 'TuftedCellsPreMuscimol.mat'])
tempPre=CumulativeCube(:,:,:,1:4);

CUMCUBEPre = cat(2,CUMCUBEPre,tempPre);
load([dataPath 'MitralCellsMuscimol.mat'])
CUMCUBEMus=CumulativeCube(:,:,:,1:4);
load([dataPath 'TuftedCellsMuscimol.mat'])
tempMus=CumulativeCube(:,:,:,1:4);
CUMCUBEMus = cat(2,CUMCUBEMus,tempMus);

%-----------------------------------------------------
ncells = size(CUMCUBEPre,2);
CellNum = 10:5:ncells;

nboot = 10;
GENERAL_PERF_Pre = zeros(20,length(CellNum),4,nboot,5);
GENERAL_PERF_Mus = zeros(20,length(CellNum),4,nboot,5);

for t = 3:20
    t
    RESPCUBE_Pre = squeeze(CUMCUBEPre(t,:,:,1:4));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 4]); 
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 4]);
    
    RESPCUBE_Mus = squeeze(CUMCUBEMus(t,:,:,1:4));
    x_mus = reshape(RESPCUBE_Mus,[ncells 5 4 4]); 
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

            CellNum(CellIter)
        end

    end
end

%%
load('TC_SVC_Poly_noRetraining')
TCperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));
TCperf_mus = squeeze(nanmean(GENERAL_PERF_Mus,3));

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

load('MC_SVC_Poly_noRetraining')
MCperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));
MCperf_mus = squeeze(nanmean(GENERAL_PERF_Mus,3));

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
%%
CellIdx = 39; % corresponds to 200 cells
dprime = [];
for odor = 1:5 
    for t = 3:20
        temp1 = squeeze(TCperf_pre(t,CellIdx,:,odor));
        temp2 = squeeze(TCperf_mus(t,CellIdx,:,odor));
        temp3 = squeeze(MCperf_pre(t,CellIdx,:,odor));
        temp4 = squeeze(MCperf_mus(t,CellIdx,:,odor));
        
        dprime(odor,t,1) = abs(nanmean(temp1)-nanmean(temp3))/sqrt(0.5*(var(temp1)+var(temp3)));
        dprime(odor,t,2) = abs(nanmean(temp1)-nanmean(temp2))/sqrt(0.5*(var(temp1)+var(temp2)));
        dprime(odor,t,3) = abs(nanmean(temp3)-nanmean(temp4))/sqrt(0.5*(var(temp3)+var(temp4)));
    end    
end
figure;hold on;
dprimeExample = squeeze(dprime(4,:,:));
plot(1:20,dprimeExample(:,1),'-k')
plot(1:20,dprimeExample(:,2),'--r')
plot(1:20,dprimeExample(:,3),'--b')
set(gca,'TickDir','out')

figure;hold on;
temp1 = squeeze(dprime(:,:,1));
cdfplot(temp1(:));
temp2 = squeeze(dprime(:,:,2));
cdfplot(temp2(:));
temp3 = squeeze(dprime(:,:,3));
cdfplot(temp3(:));
set(gca,'TickDir','out')

[h,p] = kstest(temp1(:))
[h,p] = kstest2(temp2(:),temp3(:))

figure;hold on;
shadedErrorBar(1:20,mean(temp1,1),std(temp1,0,1)./sqrt(5),'LineProps','-k')
shadedErrorBar(1:20,mean(temp2,1),std(temp2,0,1)./sqrt(5),'LineProps','--r')
shadedErrorBar(1:20,mean(temp3,1),std(temp3,0,1)./sqrt(5),'LineProps','--b')
set(gca,'TickDir','out')

figure;plot([1 2],[mean(temp2,2) mean(temp3,2)],'-o')
hold on;
errorbar([ones(5,1) 2*ones(5,1)],[mean(temp2,2) mean(temp3,2)],...
    [std(temp2,0,2)/sqrt(20) std(temp3,0,2)/sqrt(20)])
%%
figure;hold on;
CellIdx = 39;
temp1 = squeeze(TCperf_pre(:,CellIdx,:,4));
temp2 = squeeze(TCperf_mus(:,CellIdx,:,4));
temp3 = squeeze(MCperf_pre(:,CellIdx,:,4));
temp4 = squeeze(MCperf_mus(:,CellIdx,:,4));

shadedErrorBar(1:20,mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')
shadedErrorBar(1:20,mean(temp2,2),std(temp2,0,2)./sqrt(10),'LineProps','--r')
shadedErrorBar(1:20,mean(temp3,2),std(temp3,0,2)./sqrt(10),'LineProps','-b')
shadedErrorBar(1:20,mean(temp4,2),std(temp4,0,2)./sqrt(10),'LineProps','--b')
set(gca,'TickDir','out')

figure;hold on;
TimePoint = 8;
temp1 = squeeze(TCperf_pre(TimePoint,:,:,4));
temp2 = squeeze(TCperf_mus(TimePoint,:,:,4));
temp3 = squeeze(MCperf_pre(TimePoint,:,:,4));
temp4 = squeeze(MCperf_mus(TimePoint,:,:,4));

shadedErrorBar(1:size(temp1,1),mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')
shadedErrorBar(1:size(temp2,1),mean(temp2,2),std(temp2,0,2)./sqrt(10),'LineProps','--r')
shadedErrorBar(1:size(temp3,1),mean(temp3,2),std(temp3,0,2)./sqrt(10),'LineProps','-b')
shadedErrorBar(1:size(temp4,1),mean(temp4,2),std(temp4,0,2)./sqrt(10),'LineProps','--b')
set(gca,'TickDir','out')

%% Example odor 4 - dil 3
EgOdor = 4;EgConc = 3;CellIdx = 39;TimePoint = 8;

load('TC_SVC_Poly_noRetraining')
xx = squeeze(GENERAL_PERF_Pre(:,:,EgConc,:,EgOdor));
load('MC_SVC_Poly_noRetraining')
yy = squeeze(GENERAL_PERF_Pre(:,:,EgConc,:,EgOdor));

temp1 = squeeze(xx(:,CellIdx,:));
temp2 = squeeze(yy(:,CellIdx,:));
figure;
shadedErrorBar(1:size(temp1,1),mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')
shadedErrorBar(1:size(temp2,1),mean(temp2,2),std(temp2,0,2)./sqrt(10),'LineProps','-b')

temp3 = squeeze(xx(TimePoint,:,:));
temp4 = squeeze(yy(TimePoint,:,:));
figure;
shadedErrorBar(1:size(temp3,1),mean(temp3,2),std(temp3,0,2)./sqrt(10),'LineProps','-r')
shadedErrorBar(1:size(temp4,1),mean(temp4,2),std(temp4,0,2)./sqrt(10),'LineProps','-b')

dprime = [];
for t = 3:20
    temp1 = squeeze(xx(t,CellIdx,:));
    temp2 = squeeze(yy(t,CellIdx,:));   
    dprime(t) = abs(nanmean(temp1)-nanmean(temp2))/sqrt(0.5*(var(temp1)+var(temp2)));   
end    

figure;hold on;
plot(1:20,dprime,'-k')
plot(1:20,dprimeExample(:,2),'--r')
plot(1:20,dprimeExample(:,3),'--b')
set(gca,'TickDir','out')

%% Modulation Index
CellIdx = 39; % corresponds to 200 cells
modIndex = [];
for odor = 1:5 
    for t = 3:20
        temp1 = squeeze(TCperf_pre(t,CellIdx,:,odor));
        temp2 = squeeze(TCperf_mus(t,CellIdx,:,odor));
        temp3 = squeeze(MCperf_pre(t,CellIdx,:,odor));
        temp4 = squeeze(MCperf_mus(t,CellIdx,:,odor));
        
        modIndex(odor,t,1) = (nanmean(temp1)-nanmean(temp3))./(nanmean(temp1)+nanmean(temp3));
        modIndex(odor,t,2) = (nanmean(temp1)-nanmean(temp2))./(nanmean(temp1)+nanmean(temp2));
        modIndex(odor,t,3) = (nanmean(temp3)-nanmean(temp4))./(nanmean(temp3)+nanmean(temp4));
    end    
end
figure;hold on;
modIndexExample = squeeze(modIndex(4,:,:));
plot(1:20,modIndexExample(:,1),'-k')
plot(1:20,modIndexExample(:,2),'--r')
plot(1:20,modIndexExample(:,3),'--b')
set(gca,'TickDir','out')

figure;hold on;
temp1 = squeeze(modIndex(:,:,1));
cdfplot(temp1(:));
temp2 = squeeze(modIndex(:,:,2));
cdfplot(temp2(:));
temp3 = squeeze(modIndex(:,:,3));
cdfplot(temp3(:));
set(gca,'TickDir','out')
%
figure;hold on;range = 0:1:10;
temp1 = squeeze(modIndex(:,:,1));
h1 = hist(temp1(:),range);
plot(range,h1./sum(h1),'-k');
temp2 = squeeze(modIndex(:,:,2));
h2 = hist(temp2(:),range);
plot(range,h2./sum(h2),'-r');
temp3 = squeeze(modIndex(:,:,3));
h3 = hist(temp3(:),range);
plot(range,h3./sum(h3),'-b');
set(gca,'TickDir','out')

[h,p] = kstest(temp1(:))
[h,p] = kstest2(temp2(:),temp3(:))

figure;hold on;
shadedErrorBar(1:20,mean(temp1,1),std(temp1,0,1)./sqrt(5),'LineProps','-k')
shadedErrorBar(1:20,mean(temp2,1),std(temp2,0,1)./sqrt(5),'LineProps','--r')
shadedErrorBar(1:20,mean(temp3,1),std(temp3,0,1)./sqrt(5),'LineProps','--b')
set(gca,'TickDir','out')