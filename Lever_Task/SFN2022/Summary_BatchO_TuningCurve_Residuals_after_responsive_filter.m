% Summary_BatchO_TuningCurve_Residuals
% written by MD

%% USER INPUT
filter = 0; % 0 for no responsiveness filter 

%%
mousename = {'O1', 'O2', 'O3', 'O7', 'O8', 'O9'}; % dropping off O5 since noisy ephys during both open loop sessions

SessionPath = {'O1/O1_20211012_r0_processed.mat', 'O2/O2_20211011_r0_processed.mat',...
    'O3/O3_20211005_r0_processed.mat', 'O7/O7_20220630_r0_processed.mat', ...
    'O8/O8_20220702_r0_processed.mat','O9/O9_20220630_r0_processed.mat'};



%%
XVar = (0:0.0001:0.35)';
counts = zeros(length(mousename),3,4); % mouse x cutoff x comparison (CL/CL; CL/AR; CL/PR; AR/PR)
UnitsNumber = zeros(1,length(mousename));
percentmodulated = zeros(length(mousename),3,4); % mouse x cutoff x comparison (CL/CL; CL/AR; CL/PR; AR/PR)

for mouse = 1:length(mousename)
    % find responsive cells 
    switch filter
        case 1
        ThisComputerPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';
        [resptest_response_class, channels_perc] = ClassifyOdorResponsiveCells (fullfile(ThisComputerPath,SessionPath{mouse}), 1);
        cell_to_keep = find(resptest_response_class~=0);
        %get odor tuning curves
        %with version 1
        %         [TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
        %         GetOdorTuningCurves(SessionPath{mouse},cell_to_keep);
        %with version 2
        [CR_boot,TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
            GetOdorTuningCurves_v2(SessionPath{mouse},cell_to_keep);
        clearvars cell_to_keep
        case 0 
        %get odor tuning curves
        %with version 1
        %         [TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
        %         GetOdorTuningCurves(SessionPath{mouse},cell_to_keep);
        %with version 2
        [CR_boot,TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
            GetOdorTuningCurves_v2(SessionPath{mouse});
    end
    

    
    AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
    [ResidualDist] = CumDist(AllResiduals,XVar);
    cutoff(1,1) = XVar(find(ResidualDist(:,4)>=0.95,1,'first'));
    cutoff(1,2) = XVar(find(ResidualDist(:,4)>=0.97,1,'first'));
    cutoff(1,3) = XVar(find(ResidualDist(:,4)>=0.99,1,'first'));
    
    UnitsNumber(mouse) = numel(unique(PairedResiduals(:,4))); 
    
    for k = 1:3 % for each cutoff 95, 97 and 99 
        counts(mouse,k,1) = numel(unique(PairedResiduals(find(AllResiduals(:,4)>cutoff(1,k)),4)));
        percentmodulated (mouse,k,1) =  counts(mouse,k,1)*100/ UnitsNumber(mouse);
        for i = 1:3 %for each comparison 
            counts(mouse,k,i+1) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff(1,k)),4)));
            percentmodulated (mouse,k,i+1) =  counts(mouse,k,i+1)*100/ UnitsNumber(mouse);
        end
    end
end

%% plotting 
%%
figure() %number of units per mouse
X = categorical(mousename);
X = reordercats(X,mousename);
Y1 = UnitsNumber;
bar(X,Y1, 'FaceColor', Plot_Colors('k'))
xlabel('mouse id')
ylabel('# units') 
set(gca,'TickDir','out');
box off
%%
figure()%modulated 97% cutoff per mouse

Y2 = [];
chosen_cutoff = 2;
colors = [Plot_Colors('k');Plot_Colors('r');Plot_Colors('t'); Plot_Colors('o')];
for whatmouse = 1:length(X)
    Y2 = squeeze(percentmodulated(whatmouse,chosen_cutoff,:));
    b = bar(X(whatmouse),Y2); hold on
    for c = 1:length(Y2)
        b(c).FaceColor = colors(c,:);
    end 
end
hold off; 
xlabel('mouse id')
ylabel('% unit modulated - 97%cutoff')
set(gca,'TickDir','out');
lgd = legend('CL-CL','CL-AR','CL-PR', 'AR-PR');
lgd.Location = 'northeastoutside';
box off



%for whatcomp = 2:4 %CL/AR; CL/PR; AR/PR
    %Y2 = [percentmodulated(:,1,2);percentmodulated(:,1,3);percentmodulated(:,1,4)]
    
%end
%%
figure() %CL/PR - modulated with respect to number of units

scatter(Y1, percentmodulated (:,chosen_cutoff,3), 'o', 'k','filled','LineWidth',2)
refline
xlabel('# of units')
ylabel('% unit modulated during passive replay')  
box off
