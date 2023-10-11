% Summary_BatchO_Odor_Responsiveness 

mousename = {'O1', 'O2', 'O3', 'O7', 'O8', 'O9'}; % dropping off O5 since noisy ephys during both open loop sessions

MySession = {'O1/O1_20211012_r0_processed.mat', 'O2/O2_20211011_r0_processed.mat',...
    'O3/O3_20211005_r0_processed.mat', 'O7/O7_20220630_r0_processed.mat', ...
    'O8/O8_20220702_r0_processed.mat','O9/O9_20220630_r0_processed.mat'};

ThisComputerPath = 'C:\Users\Marie\Documents\data\Smellocator\Processed\Behavior';

channels_perc_batch = NaN(length(mousename),5);

for mouse = 1:length(mousename)
    SessionPath = fullfile(ThisComputerPath,MySession{mouse});
    
    [channels_perc] = ClassifyOdorResponsiveCells (SessionPath, 1);
    channels_perc_batch (mouse, :) = channels_perc; %col: u r e s m
    
end

%% plotting

figure()
X = categorical(mousename);
X = reordercats(X,mousename);
Y2 = [];
colors = [Plot_Colors('k');Plot_Colors('r');Plot_Colors('t'); Plot_Colors('o')];
for whatmouse = 1:length(X)
    Y2 = squeeze(channels_perc_batch(whatmouse,1:2));
    b = bar(X(whatmouse),Y2); hold on
    for c = 1:length(Y2)
        b(c).FaceColor = colors(c,:);
    end 
end
hold off; 
xlabel('mouse id')
ylabel('Units (%)')
ylim([20 80])
title('Percentage of responsive units');
set(gca,'TickDir','out');
lgd = legend('unresponsive','responsive');
lgd.Location = 'northeastoutside';
box off
