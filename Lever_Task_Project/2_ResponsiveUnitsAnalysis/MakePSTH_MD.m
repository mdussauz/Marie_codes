function [myFR, myPSTH, myRaster] = MakePSTH_MD (Spiketimes, Offset, WindowStart, window, varargin)
%adapted from PG MakePSTH_v3
%mostly to solve the issue of getting PSTHs of different sizes

%% parse input arguments
narginchk(1,inf)
params = inputParser;
params.CaseSensitive = false;
params.addParameter('kernelsize', 100, @(x) isnumeric(x)); %PG uses 100ms for UnitViewer
params.addParameter('downsample', 1000, @(x) isnumeric(x));


% extract values from the inputParser
params.parse(varargin{:});
kernelsize = params.Results.kernelsize;
downsample = params.Results.downsample;

Nbtrials = size(Spiketimes,1);

% Initialize raster
timeBins = window(1):window(2);
myRaster = zeros(Nbtrials,numel(timeBins)); %trials x window

if Nbtrials > 1 % if multiple trials
    for i = 1:size(Spiketimes,1)
    	% align spiketimes to the specified event
    	thisTrialSpikes = Spiketimes{i}{1} - Offset(i);
    	% convert spike times to milliseconds and floor values
    	thisTrialSpikes = ceil(1000*thisTrialSpikes);
        % remove NaNs
        thisTrialSpikes(isnan(thisTrialSpikes)) = [];
        %only keep spikes in window
        thisTrialSpikes = thisTrialSpikes(window(1)<thisTrialSpikes & thisTrialSpikes<window(2));
        % add the starting bin value to have only positive bin indices
        thisTrialSpikes = thisTrialSpikes - WindowStart; %because WindowStart should be -ve
    	% Make raster
        [C,~,ic] = unique(thisTrialSpikes);
        bin_counts = accumarray(ic,1);
        if ~isempty(C)
            % ignore any time bins less than 1 sec before window start
            bin_counts((C<=0),:) = [];
            C(C<=0) = [];
            myRaster(i,C) = bin_counts;
        end
    end
    % Make PSTH (raw)
    myPSTH = sum(myRaster,1)/i;


elseif Nbtrials ==1 % if only 1 trial

    % align spiketimes to the specified event
    thisTrialSpikes = Spiketimes - Offset;
    % convert spike times to milliseconds and floor values
    thisTrialSpikes = ceil(1000*thisTrialSpikes);
    % remove NaNs
    thisTrialSpikes(isnan(thisTrialSpikes)) = [];
    %only keep spikes in window
    thisTrialSpikes = thisTrialSpikes(window(1)<thisTrialSpikes & thisTrialSpikes<window(2));
    % add the starting bin value to have only positive bin indices
    thisTrialSpikes = thisTrialSpikes - WindowStart; %because WindowStart should be -ve
    % Make raster
    [C,~,ic] = unique(thisTrialSpikes);
    bin_counts = accumarray(ic,1);
    if ~isempty(C)
        % ignore any time bins less than 1 sec before window start
        bin_counts((C<=0),:) = [];
        C(C<=0) = [];
        myRaster(C) = bin_counts;
    end
    % Make PSTH (raw)
    myPSTH = myRaster;

end



%%

% Smoothen PSTH
taxis = -500:500;  % make a time axis of 1000 ms
gauss_kernel = normpdf(taxis, 0, kernelsize);
gauss_kernel = gauss_kernel ./ sum(gauss_kernel);

if isempty(myPSTH)
    keyboard
end

if kernelsize > 1
    myFR = 1000*conv(myPSTH,gauss_kernel,'same'); % in Hz
else
    myFR = myPSTH;
end

if downsample ~= 1000
    taxis = (1000/downsample):(1000/downsample):numel(myFR);
    if ~isempty(myFR)
        myFR = interp1(myFR,taxis');
    else
        myFR = 0*taxis;
    end
end

end