function    [f2, f3, resptest_response_class]... 
                = Responsive_bouton_ID_20221021(...
                    dataset, reference, reference_baseline, ...
                    reference_threshold, threshold, baseline_plots, ...
                    perc_plots)

%% Input - To use this code as a function, pre-allocate these variables in 
%% your source code and define the values that work for your data

% dataset % Feed with your data, the code works with data with these dimensions: frames,channel,trials

% reference                   = [56 62; 63 69; 70 76; 77 83; 84 90; 91 97];  % I am using test windows with 7 frames. It is important to use the same window size for each window you want to evaluate.
% reference_baseline          = [28 34; 35 41; 42 48; 49 55];                % You have to divide your baseline in frame windows equivalent to the ones you want to evaluate.
% reference_threshold         = [1 99; 5 95; 10 90];                         % Here you can feed possible thresholds you want to try. To get an idea of which threshold to pick, use 'Responsive_bouton_threshold_test_20221021.m'                                   
% threshold                   = 1;                                           % Pick the desired threshold by defining the used row of reference threshold.
% baseline_plots              = 1;                                           % ON/OFF option for baseline distribution plots for each channel.                                                                   
% perc_plots                  = 1;                                           % ON/OFF option for percentage anaylisis of classification.

%% Basic variables

num_Channel     = size(dataset,2);                                         % Number of ROI/Neurons/boutons/etc...
bin             = 25;                                                      % Set binning for baseline distribution plots  

xrange1         = [-2 2]; 
xlabel1         = 'zscore'; 
ylabel1         = 'counts'; 

%% Frame window averages of baseline periods. 
% Here I divided the baseline in 4 periods with a window frame size equivalent to the ones I use for the response period 
            
baseline_A(:,:) = squeeze(mean(dataset(reference_baseline(1,1):reference_baseline(1,2),:,:),1,'omitnan')); % 1st window
baseline_B(:,:) = squeeze(mean(dataset(reference_baseline(2,1):reference_baseline(2,2),:,:),1,'omitnan')); % 2nd window
baseline_C(:,:) = squeeze(mean(dataset(reference_baseline(3,1):reference_baseline(3,2),:,:),1,'omitnan')); % 3rd window
baseline_D(:,:) = squeeze(mean(dataset(reference_baseline(4,1):reference_baseline(4,2),:,:),1,'omitnan')); % 4th window

baseline_all = [baseline_A baseline_B baseline_C baseline_D];              % These values are concatenated in one single matrix

%% Baseline reference values for each channel.
               
set_threshold                   = reference_threshold(threshold,:); 

baseline_Roi_prctl_low_single   = prctile(baseline_all,set_threshold(1,1),2);    % Set the low threshold  
baseline_Roi_prctl_high_single  = prctile(baseline_all,set_threshold(1,2),2);    % Set the high threshold 

baseline_Roi_mean               = mean(baseline_all,2,'omitnan');

baseline_Roi_median             = median(baseline_all,2,'omitnan');
     

%% Baseline fistribution plots for each ROI. (Pick 'baseline_plots')

switch baseline_plots
    case 1
        f2 = figure(2); 
            for nChannel = 1:num_Channel
                    subplot(ceil(sqrt(num_Channel)),ceil(sqrt(num_Channel)),nChannel);
                        histfit(baseline_all(nChannel,:),bin,'normal');
                            title(num2str(nChannel));
                            box off;
                            set(gca,'TickDir','out');
                            set(gca,'XLim',xrange1);
                            if nChannel == 1
                                xlabel(xlabel1);
                                ylabel(ylabel1);
                            end     
                            vline(baseline_Roi_mean(nChannel,1),'b--');
                            vline(baseline_Roi_median(nChannel,1),'r--');
                            vline(baseline_Roi_prctl_low_single(nChannel,1),'k--');
                            vline(baseline_Roi_prctl_high_single(nChannel,1),'k--');                                                             
            end
    case 0
        f2 = 0;
end
  
%% Period data extraction for responsiveness test
% Currently I am dividing the after cue period in 5 parts to test if in one
% of those I found any changes. You must pick proper frame ranges by
% modifing the 'reference' array and include or remove more dataset windows
% depending on your data.

window_size = reference(1,2)-reference(1,1)+1;

dataset_periods = NaN(window_size,size(dataset,2),size(dataset,3));

for nPeriods = 1:size(reference,1)
    dataset_periods(:,:,:,nPeriods)    = dataset(reference(nPeriods,1):reference(nPeriods,2),:,:);
    
end    

% Next I reduce the frame and trial dimensions to one by getting an average value

resptest_values(:,:)    = mean(squeeze(mean(dataset_periods,1,'omitnan')),2,'omitnan'); 

%% Identification of responsive boutons, both enhanced and suppressed - SINGLE THRESHOLD VALUES

num_Category                = size(resptest_values,2);                         
resptest_index_single       = NaN(num_Channel,num_Category);

% Identification of boutons showing enhanced or suppresed responses for each period
        
for nCategory = 1:num_Category
    for nChannel = 1:num_Channel
        if resptest_values(nChannel,nCategory) > baseline_Roi_prctl_high_single(nChannel,1)
            resptest_index_single(nChannel,nCategory) = 1;

        elseif  resptest_values(nChannel,nCategory) < baseline_Roi_prctl_low_single(nChannel,1)
            resptest_index_single(nChannel,nCategory) = 2;

        else
            resptest_index_single(nChannel,nCategory) = 0;
        end
    end
end

% Classification of boutons by their responsiveness. Here I create an array
% with num Channels x 2 columns. The first column is assignated with a 1 if
% I found at least one enhanced response in one of the periods of each 
% channel row. The second with a 2 if I found at least one suppressed 
% response in each channel row. If no enhanced or suppressed response is 
% found I assign a 0 value in column 1 or 2 respectively. Then I create a 
% new variable named 'resptest_response_class' that is fed with the sum 
% value of each channel row of 'resptest_response_single'. In that way: 
% a enhanced channel corresponds to 1 + 0 = 1; 
% a suppressed channel 0 + 2 = 2; 
% a complex/dual response channel a 1 + 2 = 3; 
% an unresponsive channel a 0 + 0 = 0.

resptest_response_single    = NaN(num_Channel,2);

for nChannel = 1:num_Channel

    if find(resptest_index_single(nChannel,:) == 1,1,'first')
        resptest_response_single(nChannel,1) = 1;
    else
        resptest_response_single(nChannel,1) = 0;
    end

    if find(resptest_index_single(nChannel,:) == 2,1,'first')
        resptest_response_single(nChannel,2) = 2;
    else
        resptest_response_single(nChannel,2) = 0;
    end
end    

resptest_response_class(:,1) = resptest_response_single(:,1) + resptest_response_single(:,2);
        
%% Percentage analysis of classified channels

switch perc_plots
    case 1
        channels_u      = sum(resptest_response_class == 0);
        channels_e      = sum(resptest_response_class == 1);
        channels_s      = sum(resptest_response_class == 2);
        channels_c      = sum(resptest_response_class == 3);
        channels_all    = num_Channel;

        channels_perc(1,1) = channels_u*100/channels_all;
        channels_perc(1,3) = channels_e*100/channels_all;
        channels_perc(1,4) = channels_s*100/channels_all;
        channels_perc(1,5) = channels_c*100/channels_all;
        channels_perc(1,2) = channels_perc(1,3)+channels_perc(1,4)+channels_perc(1,5);
        
        resp_channels      = channels_e+channels_s+channels_c;
        
        resp_channels_perc(1,1) = channels_e*100/resp_channels;
        resp_channels_perc(1,2) = channels_s*100/resp_channels;
        resp_channels_perc(1,3) = channels_c*100/resp_channels;

        f3 = figure(3);
            subplot(1,2,1);
                bar(channels_perc(1,1:2))
                    box off;
                    set(gca,'TickDir','out');
                    xticklabels({'unresponsive','responsive'});
                    ylabel('Units (%)')
                    ylim([0 100]);
                    title('Percentage of responsive units');
            subplot(1,2,2);
                bar(resp_channels_perc(1,1:3))
                    box off;
                    set(gca,'TickDir','out');
                    xticklabels({'enhanced','suppressed','complex'});
                    ylabel('Responsive units (%)')
                    ylim([0 100]);
                    title('Percentage of classification of responses');
    case 0    
        f3 = 0;
end    




    
