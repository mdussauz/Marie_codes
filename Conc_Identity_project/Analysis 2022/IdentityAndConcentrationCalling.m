function [GENERAL_PERF] = IdentityAndConcentrationCalling(smoothPSTH)
% Odor Identity and concentration calling: similar to the generalization 
% across concentrations except that the classifier neuron for a given odor 
% has to fire in proportion (log scale) to the concentration without 
% explicitly requiring invariance across concentrations of the same odorant.


%% Making cumulative PSTH over relevant time period
% We want to get classification accuracy as a function
% of integration time using an ‘expanding window’ 
% consisting of feature vectors with increasing numbers of 200 ms bins 
% starting 1s before odor onset and for up to 1s after odor offset
% This section should output a matrix with 20 timewindows in dim 1
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
%cumulativePSTH(time_ind,whichcluster,whichodor,whichloc,whichrep) = FR_mean_in_bin;
                    

%% DECODER
%%
%cumulativePSTH dimensions are: time x units x odors x locations x repeats
response_matrix = CUMCUBEPre;
size(response_matrix)


%getting relevant variables `
ncells = size (response_matrix,2);
CellNum = 10:5:ncells; %performance evaluated for different numbers of neurons in steps of five.
NlocFinal = 4; %number of locations to decode
NodorFinal = 5;
Nrep = size(response_matrix,4);

t_start = 1;
t_end = size (response_matrix,1);

Nboot = 10; %number of bootstraps 

GENERAL_PERF = zeros(t_end,length(CellNum),Nrep,Nboot); % classifier output: time points X number of neuron X repeats X bootstraps

for t = t_start:t_end %for each 200 ms bin
    disp(t)
    x = squeeze(response_matrix(t,:,:,:,:)); % units x odor x loc x repeats
    x = zscore(x(:,:),0,2); %z-scoring along odor identity to avoid overfitting of neurons with high FR
    x = reshape(x,[ncells NodorFinal NlocFinal Nrep]);
    size(x)
    
  
    for test = 1:Nrep %Cross-validated performance will be evaluated across held-out trials - here 1 repeat
        % Defining the training and test sets 
        trainset = setdiff(1:4,test); %here pick 3 out of 4 repeats as training set
        xtrain = squeeze(x(:,:,:,trainset));
        xtest = squeeze(x(:,:,:,test)); % test set is the remaining repeat
        
        [nroi,ncategory,ndil,nreptrain] = size(xtrain);
        nstim = NodorFinal*NlocFinal; %number of stimuli (number of odors x number of locations)
        nTrainrep = nreptrain; %number of repeats in training set  

        % Creating the objective matrix (desired classifier target)
        Output = zeros(NodorFinal,nstim);
        Output_cross = zeros(NodorFinal,nstim,nreptrain);
 
        for k = 1:NodorFinal
            Output(k,k:NodorFinal:end) = [1:1:NlocFinal];
        end
        Output = repmat(Output,[1 nTrainrep]); % desired classifier target
        %Dimensions are odors x number of stimuli 

        for CellIter = 1:length(CellNum) % for increasing number of cells with step of 5
            PERF = zeros(Nboot,1);            
            for boot = 1:Nboot %bootstrapping - picking cells randomly with replacement 10 different times
                CellId = randperm(ncells,CellNum(CellIter)); % new vector of cells on which train/test will be performed
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';

                for k = 1:NodorFinal % for each odor 
                    Y = Output(k,:)'; %target output for that specific odor 
                    template = templateSVM('Standardize',true); % define the model
                    % returns a support vector machine learner template 
                    % suitable for training error-correcting output code (ECOC) multiclass models
                    % here the software trains the classifier using the standardized predictor matrix
                    
                    SVMModel = fitcecoc(X,Y,'Learners',template); % model fitting
                    % an ECOC model reduces the problem of classification 
                    % with 3 or more classes to a set of binary classification problems
                    
                    % X = matrix containg the predictors (training set)
                    % Y = vector containing the labels 
 
                    % From fitting we get a set of feature weights from
                    % which we can predict or classify new data (the test
                    % set):
                    [a,~]= predict(SVMModel,Xtest(:,:));
                    label(k,:) = a;

                end
                PERF(boot) = corr2(label,Output(:,1:nstim)); %compare the prediction to the truth (target output) 
                
            end

            GENERAL_PERF(t,CellIter,test,:) = PERF(:);

        end
    end
end


%% PLOTS
temp = reshape(GENERAL_PERF,[t_end size(CellNum,2) Nrep*Nboot]);
%%
% Performance across time for all cells
figure(1)

subplot(1,2,1)
time_mean = nanmean(squeeze(temp(:,end,:)),2); %mean at each timepoint across cross-validations for all cells
time_std = nanstd(squeeze(temp(:,end,:)),0,2); %std at each timepoint across cross-validations for all cells

shadedErrorBar(1*0.2:0.2:20*0.2,time_mean*100,(time_std*100)./sqrt(Nrep*Nboot))
xlabel('Time (s)')
ylabel('Classifier performance (%)')
hold on; 
xline(1, '--', {'odor','on'},'FontSize',8)
xline(3, '--', {'odor','off'},'FontSize',8)
set(gca,'TickDir','out')

% Performance across cell number at defined timepoint 
timepoint = 10; %in sec performance at time to be plotted for increasing number of cells

subplot(1,2,2)
cell_mean = nanmean(squeeze(temp(timepoint,:,:)),2); %mean for different cell number across cross-validations at one timepoint
cell_std = nanstd(squeeze(temp(timepoint,:,:)),0,2); %std for different cell number across cross-validations at one timepoint
shadedErrorBar(1*5:5:size(CellNum,2)*5,cell_mean*100,...
    (cell_std*100)./sqrt(Nrep*Nboot))

xlabel('Number of neurons')
ylabel('Classifier performance (%)')
set(gca,'TickDir','out')


%%
figure(2)
% median performance across both number of cells and time
subplot(1,2,1)
imagesc(median(temp,3)',[0 1]) %median across cross-validations
axis('square');colormap('jet'); 
xticks([5 10 15 20])
xticklabels({'1','2','3','4'})
yticks([2 4 6 8 10])
yticklabels({'10','20','30','40','50'})
set(gca,'TickDir','out')

xlabel('Time(s)')
ylabel('Number of neurons')
xline(5, 'w--',{'odor','on'},'FontSize',8)
xline(15, 'w--', {'odor','off'},'FontSize',8)
c = colorbar; c.Label.String = 'Classifier performance';


%mean performance across both number of cells and time
subplot(1,2,2)
xx = squeeze(mean(GENERAL_PERF,4));
imagesc(mean(xx,3)',[0 1])
axis('square');colormap('jet')
xticks([5 10 15 20])
xticklabels({'1','2','3','4'})
yticks([2 4 6 8 10])
yticklabels({'10','20','30','40','50'})
set(gca,'TickDir','out')

xlabel('Time(s)')
ylabel('Number of neurons')
xline(5, 'w--',{'odor','on'},'FontSize',8)
xline(15, 'w--', {'odor','off'},'FontSize',8)
c = colorbar; c.Label.String = 'Classifier performance';
end 