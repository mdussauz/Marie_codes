function [] = PlotdPCA(smoothPSTH, W, explVar)


%% This is to plot the dPCs directly
firingRates =  smoothPSTH;
firingRatesAverage = nanmean(firingRates,5);
X = firingRatesAverage(:,:);

T = size(smoothPSTH,4);
%ncomp = size(smoothPSTH,2) * size(smoothPSTH,3);
ncomp=20;
ToSmooth = 0;
sp = 25;

dpc = W(:,1:ncomp)'*X;
x = reshape(dpc,ncomp, size(smoothPSTH,2), size(smoothPSTH,3), T);

% Change these numbers to the top 3 components you need to plot;
% this will depend on whether you are plotting identity or concentration
% subspaces

%which_comps = [3 6 7]; % for AON identity with oil
which_comps = [2 4 7]; % for AON identity without oil
which_comps = [4 7 14]; % for APC identity without oil
which_comps = [2 4 8]; % for AON E3 identity without oil


figure(103);subplot(122);hold on;axis('square');
color = [165,0,38;253,174,97;116,173,209;69,117,180;49,54,149]/255;
lw = ([0.5 1 2 4]);
a = lw;
view(-25,10) % for turning the viewing angle of 3D plots

test = 6000:8000; % odor window only

for odor = 1:size(smoothPSTH, 2)
        hold on;
        for dil = 1:size(smoothPSTH, 3)
            temp2 = squeeze(x(which_comps,odor,dil,test));
            %temp2 = squeeze(x(which_comps,odor,dil,:));
            if (ToSmooth == 1)
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
 %                set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
 %              set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
                set(h1,'Color',color(odor,:),'LineWidth',a(dil));
            else
                if size(which_comps,2)==3
                    h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
%                 set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
%                 set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                elseif size(which_comps,2)==2
                   %h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)));
                   h1 = plot(squeeze(temp2(1,:)),squeeze(temp2(2,:))); %
 %                   if only 2 PC to plot
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                end 
            end
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            if size(which_comps,2)==3
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
            end 
        end
end

figure(104);subplot(122);hold on;axis('square');
%which_comps = [4 10]; % for AON concentration with oil
% which_comps = [5 8 20]; % for AON concentration without oil
% which_comps = [2 9 18]; % for APC concentration without oil
which_comps = [5 10]; % for AON E3 concentration without oil

for odor = 1:size(smoothPSTH, 2)
        hold on;
        for dil = 1:size(smoothPSTH, 3)
            temp2 = squeeze(x(which_comps,odor,dil,test));
            %temp2 = squeeze(x(which_comps,odor,dil,:));
            if (ToSmooth == 1)
                h1 = plot3(smooth(squeeze(temp2(1,:)),sp),smooth(squeeze(temp2(2,:)),sp),smooth(squeeze(temp2(3,:)),sp));
 %                set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
 %              set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
                set(h1,'Color',color(odor,:),'LineWidth',a(dil));
            else
                if size(which_comps,2)==3
                    h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)),squeeze(temp2(3,:)));
%                 set(h1,'Color',color(odor),'LineWidth',2*lw(dil),'Marker','.','MarkerSize',10*dil-5);
%                 set(h1,'Color',color(odor,:)./a(dil),'LineWidth',2);
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                elseif size(which_comps,2)==2
                   %h1 = plot3(squeeze(temp2(1,:)),squeeze(temp2(2,:)));
                   h1 = plot(squeeze(temp2(1,:)),squeeze(temp2(2,:))); %
 %                   if only 2 PC to plot
                    set(h1,'Color',color(odor,:), 'LineStyle', '-','LineWidth',a(dil))
                end 
            end
            xlabel(['dPC' num2str(which_comps(1)) '(' num2str(explVar.componentVar(which_comps(1))) ' %)']);
            ylabel(['dPC' num2str(which_comps(2)) '(' num2str(explVar.componentVar(which_comps(2))) ' %)']);
            if size(which_comps,2)==3
            zlabel(['dPC' num2str(which_comps(3)) '(' num2str(explVar.componentVar(which_comps(3))) ' %)']);
            end 
        end
end

%% This is to plot the dPCs directly

figure(105)
DPC = reshape(dpc,[20 size(smoothPSTH,2) size(smoothPSTH,3) T]); %DPC = reshape(dpc,[20 5 4 T]);
color = 'bgcrk';
which_comps = [3 6 7];

subplot(2,1,2)
for odor = 1:size(smoothPSTH,2)
    for dil = 1:size(smoothPSTH,3)
        temp = smooth(squeeze(DPC(which_comps,odor,dil,:)),10);
        plot(temp,'Color',color(odor),'LineWidth',dil/2);
        hold on;
    end
end

figure(106)
which_comps = [4 10];

subplot(2,1,2)
for odor = 1:size(smoothPSTH,2)
    for dil = 1:size(smoothPSTH,3)
        temp = smooth(squeeze(DPC(which_comps,odor,dil,:)),10);
        plot(temp,'Color',color(odor),'LineWidth',dil/2);
        hold on;
    end
end
%% bar plot with projected variances

figure(107)
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

end 