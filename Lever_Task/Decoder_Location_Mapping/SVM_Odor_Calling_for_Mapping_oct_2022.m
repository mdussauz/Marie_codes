% SVM_Odor_Calling_for_Mapping
% Written by MD
% Based on Arka's "concentration invariant odor recognition"

% Odor classification from population neural data is performed using a 
% support vector machine based decoder with non-linear polynomial kernels. 
% The feature vectors are z-scored mean odor responses for each cell. 
% Decoding accuracy at any time point (in 200 ms bins) is evaluated based 
% on mean neuronal responses from 1s before odor onset up to that time point.

% Each classifier neuron is tasked with identifying the presence of the
% odorant (non-zero value) irrespective of the location at which it was 
% presented (-15, 0, +15)

% The cross-validated performance is evaluated on held-out trials 
% Training on 3 out of 4 rep and testing on held out rep) 

% Classification performance is calculated as a correlation of the 
% decoder output with the objective matrix


% Input: PassiveTuning.RasterOut 
% PassiveTuning is a structure created by running the LocationDataParser() function on the
% 'behavior' session and after running Preprocessing on sorted ephys session 

% RasterOut contains the spike times of each neurons for different
% conditions (0 or 1 at each timestamp)
% It is a matrix with the following dimensions: 
% Units x Odors x Locations x time x repeats 
% e.g. for E3_202255 x 4 x 13 x 25501 x 4

% Trial structure in this matrix: 
% ITI = 5s; Pre Odor = 5s; Odor = 2s; Post Odor = 3.1s; ITI= 5s

% Locations: 
% Location 1 is motor -90
% Location 6 is motor -15 
% Location 7 is motor 0
% Location 8 is motor +15 
%% USER: LOADING DATA

% load data directly if PassiveTuning structure was saved - change path
% accordingly
load('C:\Users\Marie\Documents\data\Smellocator\PassiveTuningDecoder\PassiveTuning_E3_20220713.mat')

% run LocationDataParser() if PassiveTuning wasn't saved or to test
% different mice (change behavior filename if so)
%[PassiveTuning] = LocationTuningDataParser('E3_20220713_o3.mat');

%% Raster plots for sanity check 
% uncomment and change inputs as needed 
% whichneuron = 1;
% whichlocation = 7;
% hold on; PlotLocationTuningRaster(PassiveTuning.RasterOut, whichneuron, whichunit)

%% FORMATTING THE DATA TO BE INPUTED IN DECODER
%% USER: SPECIFY PARAMETERS 

% USER - Choose to smooth or not, and time window
Tosmooth =1; % 1 to smooth PSTH and 0 to bypass smoothing
timewindow = 200; % time window over which smoothing occurs - change as needed 

% USER - Choose which locations to test around 0 (index 7)
test_loc = 1; % 1 if testing -15 and +15; 2 if testing -30, -15, +15 and +30
loc_to_keep = (7 - test_loc):1:(7 + test_loc); %do not change if you want to test around center loc

% USER - Choose to keep air stimulus or not (index 1)
odor_to_keep = [2,3,4];


%% PSTH 
PSTH = PassiveTuning.RasterOut;
[Nneurons,Nodor, Nloc,lastTS, Nrep] = size(PassiveTuning.RasterOut);

%% Smooth PSTH - to get a continuous (time-dependent) rate variable
if Tosmooth == 1

    clusterNum = 1:Nneurons;
    t_wid = timewindow;  % width of kernel (in ms)
    taxis = -(t_wid*5):(t_wid*5);  % e.g. make a time axis of 1000 ms for a window of 100 ms
    gauss_kernel = normpdf(taxis, 0, t_wid);
    gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

    for clusterIdx = 1:Nneurons % loop through each unit
        for x = 1:Nodor %for each odor identity
            for y = 1:Nloc % for each location
                for j = 1:Nrep % loop through each repeat
                    tempPSTH = squeeze(PSTH(clusterIdx,x,y,:,j)); % for one cluster, all times for one type
                    zs = conv(tempPSTH,gauss_kernel,'same'); %in ms
                    smoothPSTH(clusterIdx,x,y,:,j) = zs*1000; %converting firing rate to Hz (1 ms = 1000 Hz)
                end
            end
        end
    end

    PSTH = smoothPSTH;

end
%% Making cumulative PSTH over relevant time period
% We want to get classification accuracy as a function
% of integration time using an ‘expanding window’ 
% consisting of feature vectors with increasing numbers of 200 ms bins 
% starting 1s before odor onset and for up to 1s after odor offset

% This section should output a matrix with 20 timewindows in dim 1

for whichcluster = 1:Nneurons % loop through each neuron 
    for whichodor = 1:Nodor %for each odor identity
        for whichloc = 1: Nloc % for each location
            for whichrep = 1:Nrep % for each repeat 
                time_ind = 1;
                t_min = 9000; %starts at 9s = 1s before odor period % in pre odor period
                t_max = 13000;% ends at 13s = 1s after odor period % in post odor period
                for t = (t_min+200):200:t_max %we want to average responses over expanding window 
                    tempPSTH = squeeze(PSTH(whichcluster,whichodor,whichloc,:,whichrep)); % for one cluster, all times for one type% for one cluster, all times for one type
                    FR_mean_in_bin = mean(tempPSTH(t_min:t));
                    cumulativePSTH(time_ind,whichcluster,whichodor,whichloc,whichrep) = FR_mean_in_bin;
                    t = t+200;
                    time_ind = time_ind+1;
                end
            end
        end
        
    end
    
end

%% DECODER
%%
%cumulativePSTH dimensions are: time x units x odors x locations x repeats
response_matrix = cumulativePSTH(:,:,odor_to_keep,loc_to_keep,:); %based on user specified params

%getting relevant variables `
ncells = size (response_matrix,2);
CellNum = 10:5:ncells; %performance evaluated for different numbers of neurons in steps of five.
NlocFinal = length(loc_to_keep); %number of locations to decode
NodorFinal = length(odor_to_keep);
t_start = 1;
t_end = size (response_matrix,1);

Nboot = 10; %number of bootstraps 

GENERAL_PERF = zeros(t_end,length(CellNum),Nrep,Nboot); % classifier output: time points X number of neuron X repeats X bootstraps

for t = t_start:t_end %for a specific time
    t
    x = squeeze(response_matrix(t,:,:,:,:)); % units x odor x loc x repeats
    x = zscore(x(:,:),0,2); %z-scoring along odor identity to avoid overfitting of neurons with high FR
    x = reshape(x,[ncells NodorFinal NlocFinal Nrep]);
    
  
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
            Output(k,k:NodorFinal:end) = 1;
        end
        Output = repmat(Output,[1 nTrainrep]); % desired classifier target
        
        for CellIter = 1:length(CellNum) % for different number of cells
            PERF = zeros(Nboot,NodorFinal);            
            for boot = 1:Nboot %bootstrapping - picking randomly 10 different times
                CellId = randperm(ncells,CellNum(CellIter)); % new vector of cells on which train/test will be performed
                X = xtrain(CellId,:)';
                Xtest = xtest(CellId,:)';

                for k = 1:NodorFinal % for each odor 
                    Y = Output(k,:)'; %target output for that specific odor 
                    
                    SVMModel = fitcsvm(X,Y,'Standardize',true,'KernelFunction','Poly',...
                        'KernelScale','auto'); % model fitting
                    
                  
                    % X = matrix containg the predictors (training set)
                    % Y = vector containing the labels 
 
                    % From fitting we get a set of feature weights from
                    % which we can predict or classify new data (the test
                    % set):
                    [a,b]= predict(SVMModel,Xtest(:,:));
                    temp = reshape(a,[NodorFinal,NlocFinal]);
                    label(k,:) = mean(temp,2);

                end
                PERF(boot,1:NodorFinal) = diag(label); %change text: compare the prediction to the truth (target output) 
                GENERAL_PERF(t,CellIter,test,boot,1:NodorFinal) = squeeze(PERF(boot,:));
            end
        end
    end
    t
end


%% PLOTS
%%
% Performance across time
figure(1);
subplot(1,2,1)
xx1 = squeeze(nanmean(GENERAL_PERF,3)); %TC_pre %mean across repeats
odor = 3;
temp1 = squeeze(xx1(:,end,:,odor));% all time x all cells x all boots x odor 3

shadedErrorBar(1*0.2:0.2:20*0.2,mean(temp1,2),std(temp1,0,2)./sqrt(Nboot))
xlabel('Time (s)')
ylabel('Classifier performance (%)')
hold on; 
xline(1, '--', {'odor','on'},'FontSize',8)
xline(3, '--', {'odor','off'},'FontSize',8)
set(gca,'TickDir','out')


% Performance across cell number 
subplot(1,2,2)
timepoint = 10; %in sec performance at time to be plotted for increasing number of cells

temp1 = squeeze(xx1(timepoint,:,:,odor));
shadedErrorBar(1:size(temp1,1),mean(temp1,2),std(temp1,0,2)./sqrt(Nboot))

xlabel('Number of neurons')
ylabel('Classifier performance (%)')
set(gca,'TickDir','out')


%%
figure(2)
xx2 = squeeze(mean(GENERAL_PERF,4)); %mean across boots

figCount = 1;
for j = 1:NodorFinal
    subplot(3,2,figCount);
    imagesc(squeeze(mean(xx2(:,:,:,j),3))',[0 1])
    axis('square');
    figCount = figCount+1;
    %xlim([2.5 20.5])
    set(gca,'TickDir','out')
end
colormap('jet')
xticks([5 10 15 20])
xticklabels({'1','2','3','4'})
yticks([2 4 6 8 10])
yticklabels({'10','20','30','40','50'})

xlabel('Time(s)')
ylabel('Number of neurons')
xline(5, 'w--',{'odor','on'},'FontSize',8)
xline(15, 'w--', {'odor','off'},'FontSize',8)
c = colorbar; c.Label.String = 'Classifier performance';


