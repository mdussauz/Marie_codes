function [GENERAL_PERF] = IdentityAndConcentrationCalling(smoothPSTH)
% Odor Identity and concentration calling: similar to the generalization 
% across concentrations except that the classifier neuron for a given odor 
% has to fire in proportion (log scale) to the concentration without 
% explicitly requiring invariance across concentrations of the same odorant.


%% Getting relevant variables 
global prestim
global odorstim
global poststim

if ~exist('prestim', 'var')
    %settings for 2022 experiments: 
    prestim = 6000; 
    odorstim = 2000;
    poststim = 2000; 
end 

[ncells,nodor,nconc,timepoints,nrep] = size(smoothPSTH);
nstim = nodor * nconc; 
 
%% Computing feature vector of mean response over expanding windows in relevant time period
% We want to get classification accuracy as a function
% of integration time using an ‘expanding window’ 
% consisting of feature vectors with increasing numbers of 200 ms bins 
% starting 200ms before odor onset and for up to the end of the post-stim
% period

first_bin = -200;
last_bin = odorstim + poststim; % in ms - odor window + post-stim window

temp = reshape(smoothPSTH,[ncells,nstim,timepoints,nrep]);
bins = first_bin:200:last_bin; % 200ms window bins

nbins = length(bins)-1;

response_matrix = zeros(nbins, ncells, nstim, nrep); %pre-allocating for speed 
for k = 1:nbins
    start = prestim + bins(1); 
    stop = prestim + bins(k+1); 

    response_matrix(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
    %dimensions are: time x units x odors x locations x repeats
end
                    

%% DECODER
CellNum = 10:5:ncells; %performance evaluated for different numbers of neurons in steps of five.
nboot = 10; %number of bootstraps 
GENERAL_PERF = zeros(nbins,length(CellNum),nrep,nboot); % classifier output: time points X number of neuron X repeats X bootstraps

for t = 1:nbins %for each 200 ms bin
    disp(t)
    x = squeeze(response_matrix(t,:,:,:,:)); % units x odor x loc x repeats
    x = zscore(x(:,:),0,2); %z-scoring along odor identity to avoid overfitting of neurons with high FR
    x = reshape(x,[ncells nodor nconc nrep]);
    
  
    for test = 1:nrep %Cross-validated performance will be evaluated across held-out trials - here 1 repeat
        % Defining the training and test sets 
        trainset = setdiff(1:4,test); %here pick 3 out of 4 repeats as training set
        xtrain = squeeze(x(:,:,:,trainset));
        xtest = squeeze(x(:,:,:,test)); % test set is the remaining repeat
        
        [nroi,ncategory,ndil,nreptrain] = size(xtrain);
        nstim = nodor*nconc; %number of stimuli (number of odors x number of locations)
        nTrainrep = nreptrain; %number of repeats in training set  

        % Creating the objective matrix (desired classifier target)
        Output = zeros(nodor,nstim);
        Output_cross = zeros(nodor,nstim,nreptrain);
 
        for k = 1:nodor
            Output(k,k:nodor:end) = [1:1:nconc];
        end
        Output = repmat(Output,[1 nTrainrep]); % desired classifier target
        %Dimensions are odors x number of stimuli 

        for CellIter = 1:length(CellNum) % for increasing number of cells with step of 5
            PERF = zeros(nboot,1);            
            for boot = 1:nboot %bootstrapping - picking cells randomly with replacement 10 different times
                CellId = randperm(ncells,CellNum(CellIter)); % new vector of cells on which train/test will be performed
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';

                for k = 1:nodor % for each odor 
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
temp = reshape(GENERAL_PERF,[nbins size(CellNum,2) nrep*nboot]);
%%
% Performance across time for all cells
figure()

subplot(1,2,1)
time_mean = nanmean(squeeze(temp(:,end,:)),2); %mean at each timepoint across cross-validations for all cells
time_std = nanstd(squeeze(temp(:,end,:)),0,2); %std at each timepoint across cross-validations for all cells

shadedErrorBar(bins(1:nbins)/1000,time_mean*100,(time_std*100)./sqrt(nrep*nboot))
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
shadedErrorBar(CellNum,cell_mean*100,...
    (cell_std*100)./sqrt(nrep*nboot))

xlabel('Number of neurons')
ylabel('Classifier performance (%)')
set(gca,'TickDir','out')


%%
figure()
% median performance across both number of cells and time
subplot(1,2,1)
imagesc(median(temp,3)',[0 1]) %median across cross-validations
axis('square');colormap('jet'); 
xticks(5:5:nbins)
x_labels = 1:last_bin/1000;
xticklabels(string(x_labels))
yticks(2:2:CellNum(end))
y_labels = 10:5:CellNum(end)*5;
yticklabels(string(y_labels))
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
yticks(2:2:CellNum(end))
y_labels = 10:5:CellNum(end)*5;
yticklabels(string(y_labels))
set(gca,'TickDir','out')

xlabel('Time(s)')
ylabel('Number of neurons')
% xline(5, 'w--',{'odor','on'},'FontSize',8)
% xline(15, 'w--', {'odor','off'},'FontSize',8)
c = colorbar; c.Label.String = 'Classifier performance';
end 