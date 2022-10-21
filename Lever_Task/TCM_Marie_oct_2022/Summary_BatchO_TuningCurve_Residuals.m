% Summary_BatchO_TuningCurve_Residuals
% written by MD

mousename = {'O1', 'O2', 'O3', 'O7', 'O8', 'O9'};

SessionPath = {'O1/O1_20211012_r0_processed.mat', 'O2/O2_20211011_r0_processed.mat',...
    'O3/O3_20211005_r0_processed.mat', 'O7/O7_20220630_r0_processed.mat', ...
    'O8/O8_20220702_r0_processed.mat','O9/O9_20220630_r0_processed.mat'};

XVar = (0:0.0001:0.35)';
counts = zeros(length(mousename),3,4); % mouse x cutoff x comparison (CL/CL; CL/AR; CL/PR; AR/PR)
UnitsNumber = zeros(1,length(mousename));
percentmodulated = zeros(length(mousename),3,4); % mouse x cutoff x comparison (CL/CL; CL/AR; CL/PR; AR/PR)

for mouse = 1:length(mousename)
    [TuningCurve, XBins, PairedCorrs, PairedResiduals, ControlCorrs, ControlResiduals] = ...
        GetOdorTuningCurves(SessionPath{mouse});
    
    AllResiduals = [PairedResiduals(:,1:3) ControlResiduals(:,1)];
    [ResidualDist] = CumDist(AllResiduals,XVar);
    cutoff(1,1) = XVar(find(ResidualDist(:,4)>=0.95,1,'first'));
    cutoff(1,2) = XVar(find(ResidualDist(:,4)>=0.97,1,'first'));
    cutoff(1,3) = XVar(find(ResidualDist(:,4)>=0.99,1,'first'));
    
    UnitsNumber(mouse) = PairedResiduals(end, 4); 
    
    for k = 1:3 % for each cutoff 95, 97 and 99 
        for i = 1:4 %for each comparison 
            counts(mouse,k,i) = numel(unique(PairedResiduals(find(AllResiduals(:,i)>cutoff(1,k)),4)));
            percentmodulated (mouse,k,i) =  counts(mouse,k,i)*100/ UnitsNumber(mouse);
        end
    end
end

%% plotting 
%%
figure(1) %number of units per mouse
X = categorical(mousename);
X = reordercats(X,mousename);
Y1 = UnitsNumber;
bar(X,Y1, 'FaceColor', Plot_Colors('b'))
xlabel('mouse name')
ylabel('# units') 
set(gca,'TickDir','out');
%%
figure(2)%modulated 99% cutoff per mouse

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
xlabel('mouse name')
ylabel('% unit modulated - 97%cutoff')
set(gca,'TickDir','out');
lgd = legend('CL-CL','CL-AR','CL-PR', 'AR-PR');
lgd.Location = 'northeastoutside';


%for whatcomp = 2:4 %CL/AR; CL/PR; AR/PR
    %Y2 = [percentmodulated(:,1,2);percentmodulated(:,1,3);percentmodulated(:,1,4)]
    
%end
%%
figure(3) %CL/PR - modulated with respect to number of units

scatter(Y1, percentmodulated (:,chosen_cutoff,3), 'o', 'k')
refline
xlabel('# of units')
ylabel('% unit modulated during passive replay')


    