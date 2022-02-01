
%% --------------------------------------------------------------------------------------------------
% FOR PAPER - FINAL CODES

% GENERALIZATION TO NOVEL CONC - Without Retraining 

clear;clc;close aLL;
dataPath = '/Users/abandyop/Documents/MATLAB/PhD/Imaging/AwakeOBconc/PhotonCerber/20180927/';
load([dataPath 'TuftedCellsPreMuscimol.mat'])
permuteIdx = randperm(20);

CUMCUBEPre=CumulativeCube(:,:,permuteIdx,1:4);

load([dataPath 'TuftedCellsMuscimol.mat'])
CUMCUBEMus=CumulativeCube(:,:,permuteIdx,1:4);

Results = struct('Params',{},'OutputFit',{},'OutputCross',{}, ...
    'WeightsFit',{}, 'WeightsCross', {});

ncells = size(CumulativeCube,2);
CellNum = 10:5:ncells;

GENERAL_PERF_Pre = zeros(20,length(CellNum),4,5);
GENERAL_PERF_Mus = zeros(20,length(CellNum),4,5);

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
            end
            GENERAL_PERF_Pre(t,CellIter,test,1:5) = mean(PERF(:,1,:));
            GENERAL_PERF_Mus(t,CellIter,test,1:5) = mean(PERF(:,2,:));
            CellNum(CellIter)
        end

    end
end
%%
% load('TuftedCells_NoRetraining_permutedOdors.mat')
load('TuftedCells_NoRetraining.mat')

figure;
figCount = 1;
for j = 3
    subplot(2,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Pre(:,:,j,:),4))',[0 1])
    axis('square');
    colormap('jet');xlim([3 20])
    figCount = figCount+1;
    subplot(2,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Mus(:,:,j,:),4))',[0 1])
    axis('square');
    colormap('jet');xlim([3 20])    
    figCount = figCount+1;
    temp1 = squeeze(GENERAL_PERF_Pre(:,end,j,:));
    temp2 = squeeze(GENERAL_PERF_Mus(:,end,j,:));
end
% load('MitralCells_NoRetraining_permutedOdors.mat')

load('MitralCells_NoRetraining.mat')


for j = 3
    subplot(2,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Pre(:,:,j,:),4))',[0 1])
    axis('square');
    colormap('jet');xlim([3 20])
    figCount = figCount+1;
    subplot(2,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Mus(:,:,j,:),4))',[0 1])
    axis('square');
    colormap('jet');xlim([3 20])
    figCount = figCount+1;
    temp3 = squeeze(GENERAL_PERF_Pre(:,end,j,:));
    temp4 = squeeze(GENERAL_PERF_Mus(:,end,j,:));
end

figure;hold on;
timePoint = 8;
shadedErrorBar(1:20,mean(temp1,2),std(temp1,0,2)./sqrt(5),'-r')
shadedErrorBar(1:20,mean(temp2,2),std(temp2,0,2)./sqrt(5),'--r')
shadedErrorBar(1:20,mean(temp3,2),std(temp3,0,2)./sqrt(5),'-g')
shadedErrorBar(1:20,mean(temp4,2),std(temp4,0,2)./sqrt(5),'--g')

figure;hold on;
plot(1,mean(temp1(timePoint,:),2),'-ro');
errorbar(1,mean(temp1(timePoint,:),2),std(temp1(timePoint,:),0,2)/sqrt(5),'-r','LineWidth',2)
plot(2,mean(temp2(timePoint,:),2),'-ro');
errorbar(2,mean(temp2(timePoint,:),2),std(temp2(timePoint,:),0,2)/sqrt(5),'--r','LineWidth',2)
plot(3,mean(temp3(timePoint,:),2),'-go');
errorbar(3,mean(temp3(timePoint,:),2),std(temp3(timePoint,:),0,2)/sqrt(5),'-g','LineWidth',2)
plot(4,mean(temp4(timePoint,:),2),'-go');
errorbar(4,mean(temp4(timePoint,:),2),std(temp4(timePoint,:),0,2)/sqrt(5),'--g','LineWidth',2)

xlim([0 5]);ylim([0 1])

[h,p] = ttest(squeeze(temp1(timePoint,:)),squeeze(temp2(timePoint,:)))
[h,p] = ttest(squeeze(temp3(timePoint,:)),squeeze(temp4(timePoint,:)))
[h,p] = ttest(squeeze(temp1(timePoint,:)),squeeze(temp3(timePoint,:)))

%------------------------------------------------------------------    
%% SVM for conc invariant identification

% clear;clc;close aLL;

% dataPath = '/Users/abandyop/Documents/MATLAB/PhD/Imaging/AwakeOBconc/PhotonCerber/20180927/';
% load([dataPath 'MitralCellsPreMuscimol.mat'])
% % permuteIdx = randperm(20);
% CUMCUBEPre=CumulativeCube(:,BalancedMCIdx,:,1:4);
% 
% load([dataPath 'MitralCellsMuscimol.mat'])
% CUMCUBEMus=CumulativeCube(:,BalancedMCIdx,:,1:4);

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
% ----- plotting
figure;
figCount = 1;
for j = 1:5
    subplot(5,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Pre(:,:,:,j),3))',[0 1])
%     axis('square');
    figCount = figCount+1;
    subplot(5,2,figCount);
    imagesc(squeeze(mean(GENERAL_PERF_Mus(:,:,:,j),3))',[0 1])
%     axis('square');
    figCount = figCount+1;
end
%% ------- plotting
clc

% load('TuftedCells_NoRetraining_permutedOdors.mat')
% load('TuftedCells_NoRetraining_balancedIdx_200423.mat')

% load('TuftedCells_NoRetraining.mat')

temp1 = reshape(GENERAL_PERF_Pre,[20 CellIter 4*5]);
figure;subplot(2,2,1);imagesc(mean(temp1,3)',[0 1]);axis('square')
colormap('jet');xlim([3 20]);hold on;
perf_temp = mean(temp1,3)';
for j =1:size(perf_temp,1)
    idx = find(perf_temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
    end
end
%

temp2 = reshape(GENERAL_PERF_Mus,[20 CellIter 4*5]);
subplot(2,2,2);imagesc(mean(temp2,3)',[0 1]);axis('square');
colormap('jet');xlim([3 20]);hold on;
perf_temp = mean(temp2,3)';
for j =1:size(perf_temp,1)
    idx = find(perf_temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
    end
end

% load('MitralCells_NoRetraining_balancedIdx_200423.mat')

% load('MitralCells_NoRetraining.mat')


temp3 = reshape(GENERAL_PERF_Pre,[20 CellIter 4*5]);
subplot(2,2,3);imagesc(mean(temp3,3)',[0 1]);axis('square')
colormap('jet');xlim([3 20]);hold on;
perf_temp = mean(temp3,3)';
for j =1:size(perf_temp,1)
    idx = find(perf_temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
    end
end


temp4 = reshape(GENERAL_PERF_Mus,[20 CellIter 4*5]);
subplot(2,2,4);imagesc(mean(temp4,3)',[0 1]);axis('square');
colormap('jet');xlim([3 20]);hold on;
perf_temp = mean(temp4,3)';
for j =1:size(perf_temp,1)
    idx = find(perf_temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
    end
end

timePoint = 8;
AvgPerf(1,1) = mean(squeeze(temp1(timePoint,end,:)));
AvgPerf(1,2) = std(squeeze(temp1(timePoint,end,:)));
AvgPerf(2,1) = mean(squeeze(temp2(timePoint,end,:)));
AvgPerf(2,2) = std(squeeze(temp2(timePoint,end,:)));
AvgPerf(3,1) = mean(squeeze(temp3(timePoint,end,:)));
AvgPerf(3,2) = std(squeeze(temp3(timePoint,end,:)));
AvgPerf(4,1) = mean(squeeze(temp4(timePoint,end,:)));
AvgPerf(4,2) = std(squeeze(temp4(timePoint,end,:)));


figure;hold on;
plot(1,AvgPerf(1,1),'-ro');
errorbar(1,AvgPerf(1,1),AvgPerf(1,2)/sqrt(20),'-r','LineWidth',2)
plot(2,AvgPerf(2,1),'--ro');
errorbar(2,AvgPerf(2,1),AvgPerf(2,2)/sqrt(20),'--r','LineWidth',1)
plot(3,AvgPerf(3,1),'-go');
errorbar(3,AvgPerf(3,1),AvgPerf(3,2)/sqrt(20),'-g','LineWidth',2)
plot(4,AvgPerf(4,1),'--go');
errorbar(4,AvgPerf(4,1),AvgPerf(4,2)/sqrt(20),'--g','LineWidth',1)
xlim([0 5]);ylim([0 1])

[h,p] = ttest(squeeze(temp1(timePoint,end,:,:)),squeeze(temp2(timePoint,end,:,:)))
[h,p] = ttest(squeeze(temp3(timePoint,end,:,:)),squeeze(temp4(timePoint,end,:,:)))
[h,p] = ttest2(squeeze(temp1(timePoint,end,:,:)),squeeze(temp3(timePoint,end,:,:)))


figure;hold on;
shadedErrorBar(1:20,mean(squeeze(temp1(:,end,:)),2),std(squeeze(temp1(:,end,:)),0,2)./sqrt(20),'lineProps','-r')
% shadedErrorBar(1:20,mean(squeeze(temp2(:,end,:)),2),std(squeeze(temp2(:,end,:)),0,2)./sqrt(20),'lineProps','--r')
shadedErrorBar(1:20,mean(squeeze(temp3(:,end,:)),2),std(squeeze(temp3(:,end,:)),0,2)./sqrt(20),'lineProps','-g')
% shadedErrorBar(1:20,mean(squeeze(temp4(:,end,:)),2),std(squeeze(temp4(:,end,:)),0,2)./sqrt(20),'lineProps','--g')

figure;hold on;
shadedErrorBar(1:90,mean(squeeze(temp1(end,:,:)),2),std(squeeze(temp1(end,:,:)),0,2)./sqrt(20),'lineProps','-r')
% shadedErrorBar(1:90,mean(squeeze(temp2(end,:,:)),2),std(squeeze(temp2(end,:,:)),0,2)./sqrt(20),'lineProps','--r')
shadedErrorBar(1:88,mean(squeeze(temp3(end,:,:)),2),std(squeeze(temp3(end,:,:)),0,2)./sqrt(20),'lineProps','-g')
% shadedErrorBar(1:88,mean(squeeze(temp4(end,:,:)),2),std(squeeze(temp4(end,:,:)),0,2)./sqrt(20),'lineProps','--g')

%%
temp1 = squeeze(GENERAL_PERF_Pre(:,end,:,5));
