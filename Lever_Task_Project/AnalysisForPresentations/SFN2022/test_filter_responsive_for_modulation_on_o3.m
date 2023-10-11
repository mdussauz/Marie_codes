% Test O3 odor responsive only cl ol modulation

%% USER - select which mouse open loop session to analyze
mousename = 'O3';
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
SessionPath = 'O3/O3_20211005_r0_processed.mat';
ChosenUnits = []; %ChosenUnits = [8 21 28 55 39]; %MyUnits = [8 35 28 55 39]; 

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

%% find responsive cells in closed loop
ThisComputerPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
[resptest_response_class, channels_perc] = ClassifyOdorResponsiveCells (fullfile(ThisComputerPath,SessionPath), 1);
cl_to_keep = find(resptest_response_class~=0);
cell_to_keep = cl_to_keep;
%% find responsive cells in open loop - not working right now
% [resptest_response_class,channels_perc] = ClassifyOdorResponsiveCellsReplay(fullfile(ThisComputerPath,SessionPath), 1);
% ol_to_keep = find(resptest_response_class~=0);
%cell_to_keep = unique([cl_to_keep; ol_to_keep]);
%%

%% get tuning curves
% set binsize in SmellocatorTuning.m
% use fine bins (10) or coarse bins (24)
% fine bins give better tuning curves, but noisy residual comparisons

% version 1 without bootstrap 
% [TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
%     GetOdorTuningCurves(SessionPath,cell_to_keep);

% version 2 with bootstrap 
[CR_boot,TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
    GetOdorTuningCurves_v2(SessionPath,cell_to_keep);



%% get cumulative distributions of the pairwise tuning curve correlations
XVar = (0:0.0001:0.35)';
% PairedResiduals contains in col 1 to 3 the following comparisons:
% CL-AR ; CL-PR ; AR-PR
% contains the unit number in column 4
% contains the odor number in column 5
% ControlResiduals contains the following comparison:
% CL-CL (half session vs other half picked randomly)

AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
%AllResidualsToKeep = AllResiduals(AllResiduals(:,4) == cell_to_keep);  % verify if correct
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

%% figure for SFN

figure
hold on
axis square
plot(XVar,ResidualDist(:,4),'Color',Plot_Colors('k'),'Linewidth',1.5); % CL vs CL
plot(XVar,ResidualDist(:,1),'Color',Plot_Colors('t'),'Linewidth',1.5); % CL vs AR
plot(XVar,ResidualDist(:,2),'Color',Plot_Colors('r'),'Linewidth',1.5); % CL vs PR
plot(XVar,ResidualDist(:,3),'Color',Plot_Colors('r'),'Linewidth',1.5); %AR vs PR
% line([0 0.35],[0.95 0.95],'Linestyle',':','Color','k');
% line(cutoff*[1 1], [0 1],'Linestyle',':','Color','k');

% Modulated units @95% cutoff
% for i = 1:4 %for each comparison for cutoff 95% (counts 1)
%     counts(1,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff),4)));
% end

% Modulated units @97% cutoff for cutoff 97% (counts 2)
%cutoff97 = XVar(find(ResidualDist(:,4)>=0.97,1,'first'));
line([0 0.35],[0.97 0.97],'Linestyle',':','Color','k');
line(cutoff97*[1 1], [0 1],'Linestyle',':','Color','k');
for i = 1:4 %for each comparison
    counts(2,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff97),4)));
end

% Modulated units @99% cutoff for cutoff 99% (counts 3)
%cutoff99 = XVar(find(ResidualDist(:,4)>=0.99,1,'first'));
% line([0 0.35],[0.99 0.99],'Linestyle',':','Color','r');
% line(cutoff99*[1 1], [0 1],'Linestyle',':','Color','r');
% for i = 1:4 %for each comparison
%     counts(3,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff99),4)));
% end

%Lines for 

%title(mat2str([counts(1,:) counts(2,:) counts(3,:)]))

set(gca,'TickDir','out');
lgd = legend('CL-CL','CL-AR','CL-PR', 'AR-PR');
lgd.Location = 'northeastoutside';