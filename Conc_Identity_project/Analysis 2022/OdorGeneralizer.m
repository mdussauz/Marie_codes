function [GENERAL_PERF] = OdorGeneralizer(smoothPSTH)
% Generalization to a novel concentration: after odor identity was learned 
% from a test set of concentration, we trained the classifier to identify 
% the odorant (one of a total five possibilities) by learning on three out 
% the four concentrations, and further tested for generalization to a novel
% held-out concentration.

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
end

%% DECODER
CellNum = 10:5:ncells; % number of cells considered in the analysis is increased systematically with step of 5
nboot = 10; % for each subset of cells considered, a bootstrap strategy is run 10 times
GENERAL_PERF = zeros(nbins,length(CellNum),nconc,nboot,nodor); %3rd dim is nb of conc because we are testing on held-out conc

% Creating the OBJECTIVE MATRIX (desired classifier target):
% hypothetical classifier neurons (1 for each odor) signal the
% presence (value = 1) of their corresponding odor for all four sampled 
% concentrations and its absence (value = 0) for all other odors 
Output = zeros(nodor,nodor*(nconc-1)); %output will be used to train the 
% SVM on the training set which contains 1 less concentration (held out for 
% testing)
for k = 1:nodor
    Output(k,k:nodor:end) = 1;
end
Output = repmat(Output,[1 nrep]);

label = zeros(nodor,nodor); % pre-allocating - what odor prediction based on neuronal responses and SVM weight

% Fitting and testing SVM MODEL:
tic
for t = 1:nbins %for each 200 ms bin
    disp(t)
    RESPCUBE_Pre = squeeze(response_matrix(t,:,:,:));
    x = reshape(RESPCUBE_Pre,[ncells nodor nconc nrep]); 
    x = zscore(x(:,:),0,2); %  z-scored mean odor responses for each cell (dim2 of vector x(ncells,nodor*nconc*nrep)
    x = reshape(x,[ncells nodor nconc nrep]);
       
    for test = 1:nconc %Cross-validated performance will be evaluated across held-out trials - here 1 conc 
        % Defining the training and test sets:
        trainset = setdiff(1:nconc,test);
        xtrain = squeeze(x(:,:,trainset,:));
        xtest = squeeze(x(:,:,test,:));
             
        for CellIter = 1:length(CellNum) % for increasing number of cells with step of 5 
            PERF = zeros(nboot,nodor);            
            for boot = 1:nboot %bootstrapping - picking cells randomly with replacement 10 different times
                CellId = randperm(ncells,CellNum(CellIter)); %new vector of cells on which train/test will be performed
                
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';
                
                for k = 1:nodor
                    Y = Output(k,:)';
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto');
                    
                    % X = matrix containg the predictors (training set)
                    % Y = vector containing the labels (objective matrix)
 
                    % From fitting we get a set of feature weights from
                    % which we can predict or classify new data (the test
                    % set):
                    [a,~]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[nodor,nrep]); %hack to be able to get mean of 'correct' classification for each odor at next line
                    label(k,:) = mean(temp,2); %mean classification performance across rep for each odor 
                    

                end
                PERF(boot,1:nodor) = diag(label); %reorganize to store in performance matrix
                GENERAL_PERF(t,CellIter,test,boot,1:nodor) = squeeze(PERF(boot,:));

            end

        end

    end
end
toc
%% PLOTTING FOR DECODERS

%Mean across boots and cross-validations:
xx = squeeze(nanmean(GENERAL_PERF,[3 4])); %new dim: time x cell iteration x odor 

%% 2D decoder performance as a function of time while varying 
% number of neurons included in analysis for EACH ODOR
figure(); 
figCount = 1;
for j = 1:5
    subplot(3,2,figCount);
    imagesc(squeeze(xx(:,:,j))',[0 1]) %mean across test 
    axis('square');
    figCount = figCount+1;
    xlim([first_bin last_bin])
    set(gca,'TickDir','out')
end
colormap('jet')

%% 2D decoder performance as a function of time while varying 
% number of neurons included in analysis for all odors
figure(); 
exclude_oil = 1; % to not include oil in average performance 

if exclude_oil
    selected_odors = [1:4];
else 
    selected_odors = [1:5];
end 

imagesc(squeeze(mean(xx(:,:,selected_odors),[3 4]))',[0 1]) %mean across test and odor 1:4 (ignoring oil)
axis('square');
xlim([2.5 20.5])
set(gca,'TickDir','out')
colormap('jet')

%% Performance curves across time and across number of neurons

figure(); 
% Performance for all cells number across time
subplot(1,2,1)
time_mean = nanmean(squeeze(xx(:,end,:)),2); %mean across odors at each timepoint across cross-validations for all cells
time_std = nanstd(squeeze(xx(:,end,:)),0,2); %std across odors at each timepoint across cross-validations for all cells

shadedErrorBar(bins/1000,time_mean*100,(time_std*100)./sqrt(nrep*nboot))
xlabel('Time (s)')
ylabel('Classifier performance (%)')
hold on; 
xline(1, '--', {'odor','on'},'FontSize',8)
xline(3, '--', {'odor','off'},'FontSize',8)
set(gca,'TickDir','out')

% Performance across cell number at defined timepoint (2s)
timepoint = 10; %10 is 1.8sec - performance at time to be plotted for increasing number of cells

subplot(1,2,2)
cell_mean = nanmean(squeeze(xx(timepoint,:,:)),2); %mean across odors for different cell number across cross-validations at one timepoint
cell_std = nanstd(squeeze(xx(timepoint,:,:)),0,2); %std across odors for different cell number across cross-validations at one timepoint
shadedErrorBar(CellNum,cell_mean*100,...
    (cell_std*100)./sqrt(nrep*nboot))

xlabel('Number of neurons')
ylabel('Classifier performance (%)')
set(gca,'TickDir','out')

end 