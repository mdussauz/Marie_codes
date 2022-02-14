% RunRidgeModel

%% run full model fit for PSTHs
[~, dimBeta] = ridgeMML(spikeTrace, fullR, true); %make model fit
fullFit = fullR * dimBeta; %fit data
% stimIdx = regIdx == find(ismember(regLabels,'stim'));
% % vidFit = fullR(:, ~stimIdx) * dimBeta(~stimIdx, :); %fit data
% stimFit = fullR(:, stimIdx) * dimBeta(stimIdx, :); %fit data
% save([fPath 'fullFit.mat'], 'fullFit', 'vidFit', 'stimFit');
% save([fPath 'regData.mat'], 'fullR', 'dimBeta', 'regLabels', 'regIdx');
disp('Full model fit completed');

%% cross-validated models - create for different regressors
ridgeFolds = 10;
fullV =  crossValModel(regLabels);  %cross-validated full model.
allV = NaN(size(fullV,1), size(fullV,2), 2, size(regGroups,2), 'single');
for iRegs = 1 : size(regGroups,2)
    allV(:, :, 1, iRegs) = crossValModel(regGroups{2,iRegs}(:)'); %regressor-alone model
    allV(:, :, 2, iRegs) = crossValModel(regLabels(~ismember(regLabels,regGroups{2,iRegs}(:)'))); %regressor-exclude model
    fprintf('CrossVal for %s complete (%g/%g)\n', regGroups{1,iRegs}, iRegs, size(regGroups,2));
end
disp('Cross-validation completed');

%% plotting
for y = 1:size(spikeTrace,2)
    cFull(y,1) = corr2(spikeTrace(:,y), fullV(:,y))^2;
    for z = 1:size(allV,4)
        cRegs(y,1,z) = corr2(spikeTrace(:,y), allV(:,y,1,z))^2;
        cRegs(y,2,z) = corr2(spikeTrace(:,y), allV(:,y,2,z))^2;
    end
end

% cFull = cat(1,cFull{:});
% cRegs = cat(1,cRegs{:});
% % concatenating here might not be necessary 

OdorReg = ismember(regGroups(1,:), 'Odor');
LeverReg = ismember(regGroups(1,:), 'Lever');
MotorReg = ismember(regGroups(1,:), 'Motor');
SniffReg = ismember(regGroups(1,:), 'Sniff');
% RewardsReg = ismember(regGroups(1,:), 'Rewards');
% LicksReg = ismember(regGroups(1,:), 'Licks');

figure
subplot(1,2,1)
[~, c] = sort(cFull,'descend');
plot(cFull(c), 'linewidth', 4, 'color', 'k'); hold on;
plot(cRegs(c, 1, OdorReg), 'linewidth', 4, 'color', 'r');
plot(cRegs(c, 1, LeverReg), 'linewidth', 4, 'color', 'g');
plot(cRegs(c, 1, MotorReg), 'linewidth', 4, 'color', 'b');
xlabel('Neurons');
ylabel('expl var - cross-val. R^2');
grid on; axis square
title('Predicted R^2');
legend('full model','Odor','Lever', 'Motor')

subplot(1,2,2)
% fullMean = [mean(cRegs(:, 1, OdorReg)) mean(cRegs(:, 1, LeverReg)) mean(cRegs(:, 1, MotorReg))...
%     mean(cRegs(:, 1, SniffReg)) mean(cRegs(:, 1, RewardsReg)) mean(cRegs(:, 1, LicksReg))];
% fullError = [sem(cRegs(:, 1, OdorReg)) sem(cRegs(:, 1, LeverReg)) sem(cRegs(:, 1, MotorReg))...
%     sem(cRegs(:, 1, SniffReg)) sem(cRegs(:, 1, RewardsReg)) sem(cRegs(:, 1, LicksReg))];
fullMean = [mean(cRegs(:, 1, OdorReg)) mean(cRegs(:, 1, LeverReg)) mean(cRegs(:, 1, MotorReg)) mean(cRegs(:, 1, SniffReg))];
fullError = [sem(cRegs(:, 1, OdorReg)) sem(cRegs(:, 1, LeverReg)) sem(cRegs(:, 1, MotorReg)) sem(cRegs(:, 1, SniffReg))];

errorbar(fullMean,fullError,'k-','linestyle','none','lineWidth',3); hold on
bar(fullMean,'FaceColor','g','EdgeColor','k','BarWidth',0.5,'LineWidth',2);

% uniqueMean = [mean(cFull-cRegs(:, 2, OdorReg)) mean(cFull-cRegs(:, 2, LeverReg)) mean(cFull-cRegs(:, 2, MotorReg))...
%     mean(cFull-cRegs(:, 2, SniffReg)) mean(cFull-cRegs(:, 2, RewardsReg)) mean(cFull-cRegs(:, 2, LicksReg))];
% uniqueError = [sem(cFull-cRegs(:, 2, OdorReg)) sem(cFull-cRegs(:, 2, LeverReg)) sem(cFull-cRegs(:, 2, MotorReg))...
%     sem(cFull-cRegs(:, 2, SniffReg)) sem(cFull-cRegs(:, 2, RewardsReg)) sem(cFull-cRegs(:, 2, LicksReg))];
uniqueMean = [mean(cFull-cRegs(:, 2, OdorReg)) mean(cFull-cRegs(:, 2, LeverReg)) mean(cFull-cRegs(:, 2, MotorReg)) mean(cFull-cRegs(:, 2, SniffReg))];
uniqueError = [sem(cFull-cRegs(:, 2, OdorReg)) sem(cFull-cRegs(:, 2, LeverReg)) sem(cFull-cRegs(:, 2, MotorReg)) sem(cFull-cRegs(:, 2, SniffReg))];

errorbar(-uniqueMean,uniqueError,'k-','linestyle','none','lineWidth',3); hold on
bar(-uniqueMean,'FaceColor',[25 111 61]/255,'EdgeColor','k','BarWidth',0.5,'LineWidth',2);

ax = gca;
set(ax,'xTick',1:size(fullMean,2))
set(ax,'xTickLabel',{'Odor' 'Lever', 'Motor', 'Sniff', 'Rewards', 'Licks'})
set(ax,'XTickLabelRotation',45)
ax.TickLength = [0 0];
ylabel('cross-val. R^2 R^2'); ylim([-0.2 0.2]);
axis square

%% nested function
function [Vm, cBeta, cLabels] =  crossValModel(cLabels)
        rng(1) % for reproducibility
        Vm = zeros(size(spikeTrace),'single'); %pre-allocate reconstructed spikeTrace
        randIdx = randperm(size(spikeTrace,1)); %generate randum number index
        foldCnt = floor(size(spikeTrace,1) / ridgeFolds);
        cIdx = ismember(regIdx, find(ismember(regLabels,cLabels))); %get index for task regressors
        
        for iFolds = 1:ridgeFolds
            dataIdx = true(1,size(spikeTrace,1));
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
            if iFolds == 1
                [cRidge, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true); %get beta weights and ridge penalty for task only model
            else
                [~, cBeta] = ridgeMML(spikeTrace(dataIdx,:), fullR(dataIdx,cIdx), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
            end
            Vm(~dataIdx, :) = (fullR(~dataIdx,cIdx) * cBeta); %predict remaining data
            
            if rem(iFolds,ridgeFolds/5) == 0
                fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            end
        end
    end

end