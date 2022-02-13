function [MyPeaks, Idx, Respiration_DS_Filt] = GetRespirationTimeStamps(Respiration_DS, varargin)

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('threshold', 0.1, @(x) isnumeric(x));
params.addParameter('plotfigures', false, @(x) islogical(x) || x==0 || x==1);

% extract values from the inputParser
params.parse(varargin{:});
threshold = params.Results.threshold;
plotfigures = params.Results.plotfigures;

% filter the data
order = 2; 
width = 21;
Respiration_DS_Filt = sgolayfilt(Respiration_DS,order,width);

% % rescale the data
% Respiration_DS_Filt = Respiration_DS_Filt - median(Respiration_DS_Filt);

% detect peaks and valleys
% min distance between peaks is constrained by the sampling rate (1000 hz)
% max expected sniff freq = 20 hz = 50 samples 
[MyPeaks, Idx] = findpeaks(-Respiration_DS_Filt,'MinPeakDistance',65,'MinPeakProminence',threshold);

if plotfigures
    figure;
    subplot(2,2,[1 2]);
    plot(1:length(Respiration_DS),Respiration_DS_Filt);
    hold on
    plot(Idx,-MyPeaks,'or');
    set(gca,'XLim',[10000 15000]);
    
    
    subplot(2,2,3);
    % plot the Inter-Sniff-Interval histogram
    ISIs = diff(Idx); % in milliseconds
    H = histogram(1000./ISIs,[0:1:15 Inf]);
    myHist = H.Values(1:end-1);
    bar(myHist,1,'EdgeColor','none');
    set(gca,'Color','none','YLim',[0 100*ceil(max(myHist)/100)]);
    line([10 10],[0 100*ceil(max(myHist)/100)],'color','k');
    set(gca, 'XTick', [1:15]);
end

end