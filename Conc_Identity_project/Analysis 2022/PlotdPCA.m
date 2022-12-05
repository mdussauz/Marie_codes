function [] = PlotdPCA(smoothPSTH, W, explVar, whichMarg)
%% settings for plotting 
time_window = 6000:8000; % odor window only - time window of PCs to plot 
ToSmooth = 0;
sp = 25; % window for smoothing 

color = [165,0,38;253,174,97;116,173,209;69,117,180;49,54,149]/255; % for dPC plots
lw = ([0.5 1 2 4]); % line width for dPC plots
view(-25,10) % for turning the viewing angle of 3D plots

Nodor = size(smoothPSTH, 2); 
Nconc = size(smoothPSTH, 3);

%% Computing the variables needed for plotting 
firingRates =  smoothPSTH;
firingRatesAverage = nanmean(firingRates,5);
X = firingRatesAverage(:,:);

T = size(smoothPSTH,4); %time dimension
ncomp=20; %number of components - number fixed in RunPCA

dpc = W(:,1:ncomp)'*X; 
x = reshape(dpc,[ncomp, Nodor, Nconc, T]);

%% which_comps
% Getting the top 3 components to plot for each marginalization (id an conc);
% this will depend on whether you are plotting identity or concentration
% subspaces
% in "whichMarg": 1 is identity, 2 is concentration, 3 is concentration
% independant, and 4 is conc/id interaction


%% 3D plot - identity subspace 
which_comps = find(whichMarg == 1,3); %1st 3 id comp among the 20 components

figure();subplot(122);hold on;axis('square');

for odor = 1:Nodor
        hold on;
        for dil = 1:Nconc
            temp2 = squeeze(x(which_comps,odor,dil,time_window));
            if (ToSmooth == 1) %if smoothing dCPs
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
                set(h1,'Color',color(odor,:),'LineWidth',lw(dil));
            else
                if size(which_comps,2)==3 %if 3 PCs to plot
                    h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',lw(dil))
                elseif size(which_comps,2)==2 %if 2 PCs to plot
                   h1 = plot(squeeze(temp2(1,:)),squeeze(temp2(2,:)));            
                   set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',lw(dil))
                end 
            end
            
            %labelling axis: 
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            if size(which_comps,2)==3
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
            end 
        end
end

%% 3D plot - concentration subspace 
figure();subplot(122);hold on;axis('square');
which_comps = find(whichMarg == 2,3); %1st 3 conc comp. among the 20 components

for odor = 1:Nodor
        hold on;
        for dil = 1:Nconc
            temp2 = squeeze(x(which_comps,odor,dil,time_window));
            if (ToSmooth == 1)
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
                set(h1,'Color',color(odor,:),'LineWidth',lw(dil));
            else
                if size(which_comps,2)==3 %if 3 PCs to plot
                    h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',lw(dil))
                elseif size(which_comps,2)==2 %if 2 PCs to plot
                   h1 = plot(squeeze(temp2(1,:)),squeeze(temp2(2,:)));                    
                   set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',lw(dil))
                end 
            end
            
            %labelling each axis: 
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            if size(which_comps,2)==3
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
            end 
        end
end

%% To plot the dPCs directly - identity subspace
figure()
which_comps = find(whichMarg == 1,3); %1st three id comp. among the 20 components

subplot(2,1,2)
for odor = 1:Nodor
    for dil = 1:Nconc
        temp = smooth(squeeze(x(which_comps,odor,dil,:)),10);
        plot(temp,'Color',color(odor,:),'LineWidth',dil/2);
        hold on;
    end
end

%% To plot the dPCs directly - concentration subspace
figure()
which_comps = find(whichMarg == 2,3); %1st three conc comp. among the 20 components

subplot(2,1,2)
for odor = 1:Nodor
    for dil = 1:Nconc
        temp = smooth(squeeze(x(which_comps,odor,dil,:)),10);
        plot(temp,'Color',color(odor,:),'LineWidth',dil/2);
        hold on;
    end
end

%% bar plot with projected variances
figure()
numCompToShow =15;
margColours = [250 140 0; 0 115 255; 0 0 0; 150 150 150]/256;

if ~isempty(explVar)
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

end 