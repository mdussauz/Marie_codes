%clear;clc;close aLL;
response_matrix = cumulativePSTH;

ncells = size(response_matrix,2);
CellNum = 10:5:ncells;

nboot = 10;

GENERAL_PERF = zeros(20,length(CellNum),4,nboot);

t_min = 1;
t_max = size (response_matrix,1);
max_rep = 5;

for t = t_min:t_max
    t
    RESPCUBE_Pre = squeeze(response_matrix(t,:,:,1:max_rep));
    x = reshape(RESPCUBE_Pre,[ncells 5 4 max_rep]); 
    x = zscore(x(:,:),0,2);
    x = reshape(x,[ncells 5 4 max_rep]);
    
    for test = 1:max_rep
        trainset = setdiff(1:max_rep,test);
        xtrain = squeeze(x(:,:,:,trainset));

        xtest = squeeze(x(:,:,:,test));
        
        [nroi,ncategory,ndil,nrep] = size(xtrain);
        nstim = ncategory*ndil;
        nTrainRep = nrep;
        nboot = 10;
 
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
%             GENERAL_PERF_Pre(t,CellIter,test) = mean(PERF(:,1));
%             GENERAL_PERF_Mus(t,CellIter,test) = mean(PERF(:,2));
            GENERAL_PERF(t,CellIter,test,:) = squeeze(PERF(:,1));

        end
    end
    t
end


%%
figure;hold on;
temp1 = reshape(GENERAL_PERF,[20 size(CellNum,2) 4*10]);
%temp2 = reshape(GENERAL_PERF_Mus,[20 size(CellNum,2) 4*10]);

shadedErrorBar(1:20,nanmean(squeeze(temp1(:,39,:)),2),...
    nanstd(squeeze(temp1(:,39,:)),0,2)./sqrt(40),'LineProps','-b')
% shadedErrorBar(1:20,nanmean(squeeze(temp2(:,39,:)),2),...
%     nanstd(squeeze(temp2(:,39,:)),0,2)./sqrt(40),'LineProps','--b')

%%
figure;hold on;
temp1 = reshape(GENERAL_PERF,[20 size(CellNum,2) 4*10]);
% temp2 = reshape(GENERAL_PERF_Mus,[20 size(CellNum,2) 4*10]);

shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp1(8,:,:)),2),...
    nanstd(squeeze(temp1(8,:,:)),0,2)./sqrt(40),'LineProps','-b')
% shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp2(8,:,:)),2),...
%     nanstd(squeeze(temp2(8,:,:)),0,2)./sqrt(40),'LineProps','--b')




%%
figure;
% subplot(122);imagesc(mean(GENERAL_PERF_Mus,3)',[0 1])
% axis('square');colormap('jet')

subplot(121);imagesc(mean(GENERAL_PERF,3)',[0 1])
axis('square');colormap('jet')
%%
figure;
% subplot(211);
hold on;
plot(nanmean(squeeze(GENERAL_PERF(:,end,:)),2),'-b','LineWidth',2)
% plot(nanmean(squeeze(GENERAL_PERF_Mus(:,end,:)),2),'--b','LineWidth',2)
%%
% subplot(211);
hold on;
plot(nanmean(squeeze(GENERAL_PERF(:,end,:)),2),'-b')
% plot(nanmean(squeeze(GENERAL_PERF_Mus(:,end,:)),2),'--b')

%%
%load('TC_SVM_Poly_200527')
TCperf_pre = reshape(GENERAL_PERF,[20 size(GENERAL_PERF,2) 4*10]);
%TCperf_mus = reshape(GENERAL_PERF_Mus,[20 size(GENERAL_PERF_Mus,2) 4*10]);

% load('MC_SVM_Poly_200527')
% MCperf_pre = reshape(GENERAL_PERF,[20 size(GENERAL_PERF,2) 4*10]);
% MCperf_mus = reshape(GENERAL_PERF_Mus,[20 size(GENERAL_PERF_Mus,2) 4*10]);

%%
% CellIdx = 40; % corresponds to 200 cells
% modIndex = [];
% for t = 3:20
%         temp1 = squeeze(TCperf_pre(t,CellIdx,:));
% %         temp2 = squeeze(TCperf_mus(t,CellIdx,:));
% %         temp3 = squeeze(MCperf_pre(t,CellIdx,:));
% %         temp4 = squeeze(MCperf_mus(t,CellIdx,:));
%         
%         modIndex(t,1) = abs(nanmean(temp1)-nanmean(temp3))/sqrt(0.5*(nanvar(temp1)+nanvar(temp3)));
%         modIndex(t,2) = abs(nanmean(temp1)-nanmean(temp2))/sqrt(0.5*(nanvar(temp1)+nanvar(temp2)));
%         modIndex(t,3) = abs(nanmean(temp3)-nanmean(temp4))/sqrt(0.5*(nanvar(temp3)+nanvar(temp4)));   
% end
% 
% 
% figure;hold on;
% temp1 = squeeze(modIndex(1:end,1));
% cdfplot(temp1(:));
% temp2 = squeeze(modIndex(1:end,2));
% cdfplot(temp2(:));
% temp3 = squeeze(modIndex(1:end,3));
% cdfplot(temp3(:));
% set(gca,'TickDir','out')
% 
% figure;hold on;range = 0:1:10;
% temp1 = squeeze(modIndex(:,1));
% h1 = hist(temp1(:),range);
% plot(range,h1./sum(h1),'-k');
% temp2 = squeeze(modIndex(:,2));
% h2 = hist(temp2(:),range);
% plot(range,h2./sum(h2),'-r');
% temp3 = squeeze(modIndex(:,3));
% h3 = hist(temp3(:),range);
% plot(range,h3./sum(h3),'-b');
% set(gca,'TickDir','out')
% 
% [h,p] = kstest(temp1(:))
% [h,p] = kstest2(temp2(:),temp3(:))
% 
% figure;hold on;
% plot(1:20,temp1,'-k')
% plot(1:20,temp2,'-r')
% plot(1:20,temp3,'-b')
% set(gca,'TickDir','out')

%%
% CellIdx = 40; % corresponds to 200 cells
% modIndex = [];
% for t = 3:20
%         temp1 = squeeze(TCperf_pre(t,CellIdx,:));
%         temp2 = squeeze(TCperf_mus(t,CellIdx,:));
%         temp3 = squeeze(MCperf_pre(t,CellIdx,:));
%         temp4 = squeeze(MCperf_mus(t,CellIdx,:));
%         
%         modIndex(t,1) = (nanmean(temp1)-nanmean(temp3))./(nanmean(temp1)+nanmean(temp3));
%         modIndex(t,2) = (nanmean(temp1)-nanmean(temp2))./(nanmean(temp1)+nanmean(temp2));
%         modIndex(t,3) = (nanmean(temp3)-nanmean(temp4))./(nanmean(temp3)+nanmean(temp4));
% end
% 
% 
% figure;hold on;
% temp1 = squeeze(modIndex(1:end,1));
% cdfplot(temp1(:));
% temp2 = squeeze(modIndex(1:end,2));
% cdfplot(temp2(:));
% temp3 = squeeze(modIndex(1:end,3));
% cdfplot(temp3(:));
% set(gca,'TickDir','out')
% 
% figure;hold on;range = 0:1:10;
% temp1 = squeeze(modIndex(:,1));
% h1 = hist(temp1(:),range);
% plot(range,h1./sum(h1),'-k');
% temp2 = squeeze(modIndex(:,2));
% h2 = hist(temp2(:),range);
% plot(range,h2./sum(h2),'-r');
% temp3 = squeeze(modIndex(:,3));
% h3 = hist(temp3(:),range);
% plot(range,h3./sum(h3),'-b');
% set(gca,'TickDir','out')
% 
% [h,p] = kstest(temp1(:))
% [h,p] = kstest2(temp2(:),temp3(:))
% 
% figure;hold on;
% plot(1:20,temp1,'-k')
% plot(1:20,temp2,'-r')
% plot(1:20,temp3,'-b')
% set(gca,'TickDir','out')
