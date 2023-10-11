function    [f1, resptest_perctest]... 
                = Responsive_bouton_threshold_test_20221021(...
                    reference, reference_baseline, dataset)

%% Input - To use this code as a function, pre-allocate these variables in 
%% your source code and define the values that work for your data

% dataset % Feed with your data, the code works with data with these dimensions: frames,channel,trials
% dataset % Feed with your data, the code works with data with these dimensions: frames,channel,trials

% reference                   = [56 62; 63 69; 70 76; 77 83; 84 90; 91 97];  % I am using test windows with 7 frames. It is important to use the same window size for each window you want to evaluate.
% reference_baseline          = [28 34; 35 41; 42 48; 49 55];                % You have to divide your baseline in frame windows equivalent to the ones you want to evaluate.

%% Basic variables

num_Channel                 = size(dataset,2);                                         

%% Frame window averages of baseline periods. 
% Here I divided the baseline in 4 periods with a window frame size equivalent to the ones I use for the response period 
            
baseline_A(:,:) = squeeze(mean(dataset(reference_baseline(1,1):reference_baseline(1,2),:,:),1,'omitnan')); % 1st window
baseline_B(:,:) = squeeze(mean(dataset(reference_baseline(2,1):reference_baseline(2,2),:,:),1,'omitnan')); % 2nd window
baseline_C(:,:) = squeeze(mean(dataset(reference_baseline(3,1):reference_baseline(3,2),:,:),1,'omitnan')); % 3rd window
baseline_D(:,:) = squeeze(mean(dataset(reference_baseline(4,1):reference_baseline(4,2),:,:),1,'omitnan')); % 4th window

baseline_all = [baseline_A baseline_B baseline_C baseline_D];              % These values are concatenated in one single matrix
  
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

%% Variable threshold vs responsive-bouton analysis ('threshold_test' option)

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

%% Baseline average frame windows distribution percentiles    

num_Category                = size(resptest_values,2); 
threshold_low               = 1:1:50; 
threshold_high              = flip(51:1:100);
num_Thresholds              = size(threshold_low,2); 

resptest_index              = NaN(num_Channel,num_Category);
resptest_response           = NaN(num_Channel,3);
resptest_perctest           = NaN(num_Channel,num_Thresholds);

% The next section loops the responsiveness evaluation through the
% whole range of possible integer percentile thresholds with respect to
% the beaseline

for nThreshold = 1:num_Thresholds
    baseline_Roi_prctl_low      = prctile(baseline_all,threshold_low(1,nThreshold),2);  % Set the lower threshold to identify suppresed responses 
    baseline_Roi_prctl_high     = prctile(baseline_all,threshold_high(1,nThreshold),2); % Set the higher threshold to identify enhanced responses       

    % Identification of responsive boutons for each threshold and
    % evaluated frame window.

    for nCategory = 1:num_Category
        for nChannel = 1:num_Channel
            if resptest_values(nChannel,nCategory) > baseline_Roi_prctl_high(nChannel,1)
                resptest_index(nChannel,nCategory) = 1;
            elseif resptest_values(nChannel,nCategory) < baseline_Roi_prctl_low(nChannel,1)
                resptest_index(nChannel,nCategory) = 2;
            else
                resptest_index(nChannel,nCategory) = 0;
            end
        end
    end

    % Identification of boutons showing no, enhanced, suppresed or both
    % enhanced and suppresed responses
    for nChannel = 1:num_Channel

        if find(resptest_index(nChannel,:)==1,1,'first')
            resptest_response(nChannel,1) = 1;
        else
            resptest_response(nChannel,1) = 0;
        end

        if find(resptest_index(nChannel,:)==2,1,'first')
            resptest_response(nChannel,2) = 2;
        else
            resptest_response(nChannel,2) = 0;
        end
    end 
    
    resptest_response_class_test(:,1) = resptest_response(:,1)+resptest_response(:,2);
    resptest_perctest(:,nThreshold) = resptest_response_class_test(:,1);
end

%% Plot 

threshold_low               = 1:1:50; 
threshold_high              = 51:1:100;
threshold_high              = flip(threshold_high);
num_thresh                  = size(threshold_low,2); 

for nThresh = 1:num_thresh
    num_boutons(1,nThresh)      = sum(resptest_perctest(:,nThresh) == 1)+sum(resptest_perctest(:,nThresh) == 2)+sum(resptest_perctest(:,nThresh) == 3);
    num_enhanced(1,nThresh)     = sum(resptest_perctest(:,nThresh) == 1);
    num_suppressed(1,nThresh)   = sum(resptest_perctest(:,nThresh) == 2);
    num_both(1,nThresh)         = sum(resptest_perctest(:,nThresh) == 3);

end    

f1 = figure(1);
    subplot(1,2,1);
        hold on
        plot(threshold_high,num_boutons,'k','LineWidth',2.5);
        plot(threshold_high,num_enhanced,'b','LineWidth',2.5);
        plot(threshold_high,num_suppressed,'r','LineWidth',2.5);
        plot(threshold_high,num_both,'g','LineWidth',2.5);
            xlabel('threshold (%)');
            ylabel('Bouton number');
            title('Threshold-dependent classification of boutons')
            xlim([50 100]);
            ylim([0 num_Channel]);
            vline(95,'k--');
            vline(99,'r--');
    subplot(1,2,2);
        hold on
        plot(threshold_high,100*num_boutons/num_Channel,'k','LineWidth',2.5);
        plot(threshold_high,100*num_enhanced/num_Channel,'b','LineWidth',2.5);
        plot(threshold_high,100*num_suppressed/num_Channel,'r','LineWidth',2.5);
        plot(threshold_high,100*num_both/num_Channel,'g','LineWidth',2.5);
            xlabel('threshold (%)');
            ylabel('Bouton number (%)');
            xlim([50 100]);
            ylim([0 100]);
            title('Threshold-dependent classification of boutons - percentage')
            vline(95,'k--');
            vline(99,'r--');
            legend({'responsive','enhanced','suppressed','complex'});