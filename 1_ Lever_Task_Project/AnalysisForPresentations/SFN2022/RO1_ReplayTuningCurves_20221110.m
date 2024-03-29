% adapted by MD 

% Goals: 
% plot frequency cumulative plot of tuning curve comparisons 
% between closed loop, active and passive replay
% and plotting of tuning curves and rasters + PSTHs for each cond

%% USER - select which mouse open loop session to analyze
mousename = 'O3';
version = 2; %which version of tuning curve comparison
%% USER - turn on/off response plotting (tuning curve, raster and psth)
% doesn't turn off the cumulative frequency plot
Toplot =0; %1 for ON

%%
switch mousename
    case 'O1'
SessionPath = 'O1/O1_20211012_r0_processed.mat';
ChosenUnits = []; %MyUnits = [8 35 28 55 39]; 

    case 'O2'
SessionPath = 'O2/O2_20211011_r0_processed.mat';
ChosenUnits = []; %MyUnits = [8 35 28 55 39]; 

    case 'O3'
SessionPath = 'O3/O3_20211005_r0_processed.mat'; %units chosen for poster
ChosenUnits = []; %ChosenUnits = [8 21 28 55 39]; 
%ChosenUnits = [8 21 34 28 29 39 44 55]; 

    case 'O5'
% error for both of these sessions
% doesn't get solved by re-running preprocessing on 10/20/2022
%SessionPath = 'O5/O5_20211005_r0_processed.mat'; %error related to OdorTTLs
% or 
SessionPath = 'O5/O5_20211006_r0_processed.mat'; %error related to motor
ChosenUnits = []; %MyUnits = [8 35 28 55 39]; 

    case 'O8'
SessionPath = 'O8/O8_20220702_r0_processed.mat';
ChosenUnits = []; %ChosenUnits = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 

    case 'O9'
SessionPath = 'O9/O9_20220630_r0_processed.mat';
ChosenUnits = []; %ChosenUnits = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 

    case 'O7'
SessionPath = 'O7/O7_20220630_r0_processed.mat';
ChosenUnits = []; %ChosenUnits = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 

    case 'PCX4'
SessionPath = 'PCX4/PCX4_20210721_r0_processed.mat';
ChosenUnits = [2 38 45 55 64]; % may be 38

end

%% get tuning curves
% set binsize in SmellocatorTuning.m
% use fine bins (10) or coarse bins (24)
% fine bins give better tuning curves, but noisy residual comparisons
switch version
    case 1
% version 1 no bootstrap
[TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves(SessionPath);
    case 2
% version 2 with bootstraps
[mean_CR_all_boot,TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves_v2(SessionPath);
end

%% get cumulative distributions of the pairwise tuning curve correlations
XVar = (0:0.0001:0.35)';
% PairedResiduals contains in col 1 to 3 the following comparisons:
% CL-AR ; CL-PR ; AR-PR
% contains the unit number in column 4
% contains the odor number in column 5
% ControlResiduals contains the following comparison:
% CL-CL (half session vs other half picked randomly)

AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
[ResidualDist] = CumDist(AllResiduals,XVar);
cutoff = XVar(find(ResidualDist(:,4)>=0.95,1,'first'));

figure
hold on
axis square
plot(XVar,ResidualDist(:,4),'Color',Plot_Colors('k'),'Linewidth',2); % CL vs CL
plot(XVar,ResidualDist(:,1),'Color',Plot_Colors('t'),'Linewidth',2); % CL vs AR
plot(XVar,ResidualDist(:,2),'Color',Plot_Colors('r'),'Linewidth',2); % CL vs PR
plot(XVar,ResidualDist(:,3),':','Color',Plot_Colors('r'),'Linewidth',2); %AR vs PR
line([0 0.35],[0.95 0.95],'Linestyle',':','Color','k');
line(cutoff*[1 1], [0 1],'Linestyle',':','Color','k');

% Modulated units @95% cutoff
for i = 1:4 %for each comparison for cutoff 95% (counts 1)
    counts(1,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff),4)));
end

% Modulated units @97% cutoff for cutoff 97% (counts 2)
cutoff97 = XVar(find(ResidualDist(:,4)>=0.97,1,'first'));
line([0 0.35],[0.97 0.97],'Linestyle',':','Color','b');
line(cutoff97*[1 1], [0 1],'Linestyle',':','Color','b');
for i = 1:4 %for each comparison
    counts(2,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff97),4)));
end

% Modulated units @99% cutoff for cutoff 99% (counts 3)
cutoff99 = XVar(find(ResidualDist(:,4)>=0.99,1,'first'));
line([0 0.35],[0.99 0.99],'Linestyle',':','Color','r');
line(cutoff99*[1 1], [0 1],'Linestyle',':','Color','r');
for i = 1:4 %for each comparison
    counts(3,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff99),4)));
end

title(mat2str([counts(1,:) counts(2,:) counts(3,:)]))

set(gca,'TickDir','out');
lgd = legend('CL-CL','CL-AR','CL-PR', 'AR-PR');
lgd.Location = 'northeastoutside';

% Help to read cumulative frequency graph:
% Each point on the graph is the comparison of the tuning curve of one cell
% for one odor across 2 conditions (or odor-cell pair)

% To read the graph you take the residual number of cl-cl based on cutoff
% and from the y axis you get how many odor-cell pairs have a residual that
% is less than the cutoff (so not modulated) - so you cannot get directly
% from the graph the actual number of cells that are modulated 

% You do not get directly the number of modulated cells, hence the
% "counts" line:
% In PairedResiduals you have all possible cell-odor pair (units x 3 odors)
% for the 3 comparisons (not the control one)
% Each colomn of AllResiduals is a comparison (including control)
% and you find for a specific cutoff how many cell-odor pair are modulated
% in column 4 which will give you a unit number associated to that
% cell-odor pair
% unique() will then make sure you only count unit once - no matter if it got
% modulated for one odor or multiple odors

if Toplot == 1
    %% Plot the tuning curves
    if isempty(ChosenUnits)
        ChosenUnits = 1:length(TuningCurve.ClosedLoopFull(1,1,:,1)); % to get all units
    end
    for unit = 1:numel(ChosenUnits)  %max(PairedResiduals(:,4))
        if mod(unit,5) == 1
            figure;
            i = 0;
        end
        i = i + 1;
        for j = 1:3
            subplot(5,3,j+3*(i-1));
            hold on;
            plot(mean(XBins,2)'-125,TuningCurve.ClosedLoopFull(:,1,ChosenUnits(unit),j),...
                'color',Plot_Colors('k'),'Linewidth',2);
            plot(mean(XBins,2)'-125,TuningCurve.OpenLoop(:,1,ChosenUnits(unit),j),...
                'color',Plot_Colors('t'),'Linewidth',2);
            plot(mean(XBins,2)'-125,TuningCurve.Passive(:,1,ChosenUnits(unit),j),...
                'color',Plot_Colors('r'),'Linewidth',2);
            
            foo = intersect(find(PairedResiduals(:,4)==ChosenUnits(unit)),find(PairedResiduals(:,5)==j));
            if any(PairedResiduals(foo,1:3)>cutoff)
                title(mat2str([ChosenUnits(unit) PairedResiduals(foo,1:3)],2),'Color','r');
            else
                title(mat2str([ChosenUnits(unit) PairedResiduals(foo,1:3)],2));
            end
        end
    end
    
    %% Plot the selected units
    PlotReplayResponses(SessionPath,ChosenUnits);
    
end