%% bar plot with projected variances
numCompToShow =15;
margColours = [250 140 0; 0 115 255; 0 0 0; 150 150 150]/256;

if ~isempty(explVar)
    %axBar = subplot(4,4,9);
    hold on
    axis([0 numCompToShow+1 0 12.5])
    xlabel('De-mixed components')
    ylabel('Variance explained (%)')
    b = bar(explVar.margVar(:,1:numCompToShow)' , 'stacked', 'BarWidth', 0.75);
    
    % fix for newer Matlab versions, May 2018, thanks to Tucker Fisher
    for idx = 1:numel(b)
        b(idx).FaceColor = margColours(idx,:);
    end    
    % end fix
    
    caxis([1 length(margColours)+256])
end

% Legend will show names for each color
legend( {'Identity','Concentration','Condition Independent',...
    'I/C interaction'}) 

legend boxoff  % Hides the legend's axes (legend border and background)

% % set 4 display names for the 4 handles
% set(b, {'DisplayName'}, {'Identity','Concentration','Condition Independent',...
%     'I/C interaction'})
% set(b(1), 'visible', 'off')    % Hides the legend's axes (legend border and background)

%% pie chart
% if ~isempty(options.explainedVar)
%     axes('position', [0.205 0.47 0.1 0.1])
%     
%     if isfield(options.explainedVar, 'totalMarginalizedVar_signal')
%         d = options.explainedVar.totalMarginalizedVar_signal / options.explainedVar.totalVar_signal * 100;
%        
%         % In some rare cases the *signal* explained variances can be
%         % negative (usually around 0 though); this means that the
%         % corresponding marginalization does not carry [almost] any signal.
%         % In order to avoid confusing pie charts, we set those to zero and
%         % rescale the others to sum to 100%.
%         if ~isempty(find(d<0, 1))
%             d(d<0) = 0;
%             d = d/sum(d)*100;
%         end
%     else
%         d = options.explainedVar.totalMarginalizedVar / options.explainedVar.totalVar * 100;
%     end
%     
%     % Rounding such that the rounded values still sum to 100%. Using
%     % "largest remainder method" of allocation
%     roundedD = floor(d);
%     while sum(roundedD) < 100
%         [~, ind] = max(d-roundedD);
%         roundedD(ind) = roundedD(ind) + 1;
%     end
%     
%     if ~isempty(options.marginalizationNames)
%         for i=1:length(d)
%             margNamesPerc{i} = [options.marginalizationNames{i} ' ' num2str(roundedD(i)) '%'];
%         end
%     else
%         for i=1:length(d)
%             margNamesPerc{i} = [num2str(roundedD(i)) '%'];
%         end
%     end
%     pie(d, ones(size(d)), margNamesPerc)
%     caxis([1 length(options.marginalizationColours) + 256])
% end