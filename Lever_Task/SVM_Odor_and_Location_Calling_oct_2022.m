%SVM_Odor_and_Location_Calling 

% Written by MD

% Input: PassiveTuning.RasterOut 
% This structure was created by running the LocationDataParser on the
% 'behavior' session after running Preprocessing on sorted ephys session 

% It is a matrix with the following dimensions: 
% Units x Odors x Locations x time x repeats 
% e.g. 55 x 4 x 13 x 25501 x 4

%% FORMATTING THE DATA TO BE INPUTED IN DECODER
%% Smooth PSTH 

%% Making cumulative PSTH - try without that first 

%% DECODER
%%
CUMCUBEPre=CumulativeCube(:,:,:,1:4);
 
ncells = size(CumulativeCube,2);
CellNum = 10:5:ncells; %evaluated for different numbers of neurons in steps of five.

nboot = 10;

GENERAL_PERF_Pre = zeros(20,length(CellNum),4,nboot); %PreMuscimol classifier output: time points X number of neuron bins X repeats X bootstraps

for t = 3:20 %ignore the first two timepoints (~400 ms because of odor latency)
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
        nboot = 10;
 
        Output = zeros(ncategory,nstim);
        Output_cross = zeros(ncategory,nstim,nrep);
 
        for k = 1:ncategory
            Output(k,k:ncategory:end) = [1 2 3 4];
        end
        Output = repmat(Output,[1 nTrainRep]); % desired classifier target

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


%% Plotting

figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[20 size(CellNum,2) 4*10]);

shadedErrorBar(1:20,nanmean(squeeze(temp1(:,87,:)),2),...
    nanstd(squeeze(temp1(:,10,:)),0,2)./sqrt(40),'-b')


%%
figure;hold on;
temp1 = reshape(GENERAL_PERF_Pre,[20 size(CellNum,2) 4*10]);

shadedErrorBar(1:size(CellNum,2),nanmean(squeeze(temp1(8,:,:)),2),...
    nanstd(squeeze(temp1(8,:,:)),0,2)./sqrt(40),'-b')


%%
figure;
subplot(121);imagesc(median(temp1,3)',[0 1])
axis('square');colormap('jet')

%%
figure;
% subplot(211);
hold on;
plot(nanmean(squeeze(GENERAL_PERF_Pre(:,end,:)),2),'-b','LineWidth',2)

%%
% subplot(211);
hold on;
plot(nanmean(squeeze(GENERAL_PERF_Pre(:,end,:)),2),'-b')

%%
TCperf_pre = reshape(GENERAL_PERF_Pre,[20 size(GENERAL_PERF_Pre,2) 4*10]);

figure;
subplot(221);
temp = squeeze(mean(GENERAL_PERF_Pre,4));
imagesc(mean(temp,3)',[0 1])
axis('square');colormap('jet')
