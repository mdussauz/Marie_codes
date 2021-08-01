%----------------------------------
% 
% Results = struct('Params',{},'OutputFit',{},'OutputCross',{}, ...
%     'WeightsFit',{}, 'WeightsCross', {});
% ---------------------------------------------

% smoothPSTH5D = binnedPSTH;
% smoothPSTH5D = smoothPSTH5D(:,:,50:70,1:5);
smoothPSTH5D = cumulativePSTH;
%smoothPSTH5D = smoothPSTH5D(45:75,:,:,1:5);
%-------------------------------------------------
ncells = size (smoothPSTH5D,2);
CellNum = 10:5:ncells;
% CellNum = 200;

t_min = 1;
t_max = 30;

nboot = 10;
GENERAL_PERF_Pre = zeros(20,length(CellNum),4,nboot,5);

tic
for t = t_min:t_max %time %each bin is 200 ms
    t
    RESPCUBE_Pre = squeeze(smoothPSTH5D(t,:,:,1:5)); %pick a specific time
    %RESPCUBE_Pre = squeeze(smoothPSTH5D(:,:,t,1:5)); %pick a specific time
    x = reshape(RESPCUBE_Pre,[ncells 5 4 5]); %nber of cells id conc rep
    x = zscore(x(:,:),0,2); %zscore in dim 2 so here for each identity 
    x = reshape(x,[ncells 5 4 5]); %nber of cells id conc rep
    
    for test = 1:5 % rep
        trainset = setdiff(1:5,test); %here pick one of the rep as training set
        %    C = setdiff(A,B) for vectors A and B, returns the values in A that
        %    are not in B with no repetitions. C will be sorted.
        xtrain = squeeze(x(:,:,:,trainset));
        
        xtest = squeeze(x(:,:,:,test));
        
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
                
                for k = 1:ncategory
                    Y = Output(k,:)';
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto'); % change Poly to Linear
                    
                    [a,b]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[ncategory,ndil]);
                    label(k,:) = mean(temp,2);
                    
%                     [a,b]= predict(SVMModel,Xtest_mus(:,:));
%                     temp = reshape(a,[ncategory,ndil]);
%                     label_mus(k,:) = mean(temp,2);
                    
                end
                PERF(boot,1,1:5) = diag(label);
                GENERAL_PERF_Pre(t,CellIter,test,boot,1:5) = squeeze(PERF(boot,1,:));
            end

        end
        
    end
end

%% PLOTTING ------------------------------------
%load('TC_SVC_Poly_noRetraining')
TCperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));

xx = squeeze(mean(GENERAL_PERF_Pre,4));

figure;
figCount = 1;
for j = 1:5
    subplot(3,2,figCount);
    imagesc(squeeze(mean(xx(:,:,:,j),3))',[0 1])
    axis('square');
    figCount = figCount+1;
    xlim([0.5 30.5])
    set(gca,'TickDir','out')
end
colormap('jet')

% %load('MC_SVC_Poly_noRetraining')
% MCperf_pre = squeeze(nanmean(GENERAL_PERF_Pre,3));
% 
% xx = squeeze(mean(GENERAL_PERF_Pre,4));

% figure;
% figCount = 1;
% for j = 1:5
%     subplot(5,2,figCount);
%     imagesc(squeeze(mean(xx(:,:,:,j),3))',[0 1])
%     axis('square');
%     figCount = figCount+1;
%     xlim([0.5 30.5])
%     set(gca,'TickDir','out')
% end
% colormap('jet')
%%
CellIdx = 39; % corresponds to 200 cells
dprime = [];
for odor = 1:5
    for t = t_min:t_max
        temp1 = squeeze(TCperf_pre(t,CellIdx,:,odor));
%         temp3 = squeeze(MCperf_pre(t,CellIdx,:,odor));
%         dprime(odor,t) = abs(nanmean(temp1)-nanmean(temp3))/sqrt(0.5*(var(temp1)+var(temp3)));

    end
end
% figure;hold on;
% dprimeExample = squeeze(dprime(4,:));
% plot(1:t_max,dprimeExample(:,1),'-k')
% plot(1:t_max,dprimeExample(:,2),'--r')
% plot(1:t_max,dprimeExample(:,3),'--b')
% set(gca,'TickDir','out')

% figure;hold on;
% temp1 = dprime(:,:);
% cdfplot(temp1(:));
% set(gca,'TickDir','out')
% 
% figure;hold on;range = 0:1:10;
% temp1 = dprime(:,:);
% h1 = hist(temp1(:),range);
% plot(range,h1./sum(h1),'-k');
% set(gca,'TickDir','out')
% 
% [h,p] = kstest(temp1(:))

% 
% figure;hold on;
% shadedErrorBar(1:t_max,mean(temp1,1),std(temp1,0,1)./sqrt(5),'LineProps','-k')
% set(gca,'TickDir','out')
%%
figure;hold on;
CellIdx = 39;
temp1 = squeeze(TCperf_pre(:,CellIdx,:,4));
% temp3 = squeeze(MCperf_pre(:,CellIdx,:,4));

shadedErrorBar(1:t_max,mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')
% shadedErrorBar(1:20,mean(temp3,2),std(temp3,0,2)./sqrt(10),'LineProps','-b')
set(gca,'TickDir','out')

figure;hold on;
TimePoint = 15; %15 for AON
%TimePoint = 8;
temp1 = squeeze(TCperf_pre(TimePoint,:,:,4));
% temp3 = squeeze(MCperf_pre(TimePoint,:,:,4));

shadedErrorBar(1:size(temp1,1),mean(temp1,2),std(temp1,0,2)./sqrt(10),'LineProps','-r')
% shadedErrorBar(1:size(temp3,1),mean(temp3,2),std(temp3,0,2)./sqrt(10),'LineProps','-b')
set(gca,'TickDir','out')

%%
% CellIdx = 39; % corresponds to 200 cells
% modIndex = [];
% for odor = 1:5
%     for t = t_min:t_max
%         temp1 = squeeze(TCperf_pre(t,CellIdx,:,odor));
%         temp3 = squeeze(MCperf_pre(t,CellIdx,:,odor));      
%         modIndex(odor,t) = (nanmean(temp1)-nanmean(temp3))./(nanmean(temp1)+nanmean(temp3));
% 
%     end
% end
% figure;hold on;
% modIndexExample = squeeze(modIndex(4,:));
% plot(1:20,modIndexExample(:,1),'-k')
% plot(1:20,modIndexExample(:,2),'--r')
% plot(1:20,modIndexExample(:,3),'--b')
% set(gca,'TickDir','out')
% 
% figure;hold on;
% temp1 = modIndex(:,:);
% cdfplot(temp1(:));
% set(gca,'TickDir','out')
% %
% figure;hold on;range = 0:1:10;
% temp1 = modIndex(:,:);
% h1 = hist(temp1(:),range);
% plot(range,h1./sum(h1),'-k');
% set(gca,'TickDir','out')
% 
% [h,p] = kstest(temp1(:))
% 
% figure;hold on;
% shadedErrorBar(1:t_max,mean(temp1,1),std(temp1,0,1)./sqrt(5),'LineProps','-k')
% set(gca,'TickDir','out')
