clear;clc;close all;
dataPath = 'C:\Users\Marie\Documents\data\Data_Photometry_Diego';

load([dataPath 'APC_post.mat';])
%Data_compilation is mouse x time x odor x rep
temp=Data_Compilation(:,:,:,1:5); 
%temp=Data_Compilation(:,:,:,1:10);% I will need to try different
%combinations of rep
bins = 0:0.2:4;
Fs = 1000;

%% Here cumulating responses for different time bins - need to try with and without cumulating
% AONPre = [];
% for k = 1:length(bins)-1
%     start = Fs*(10 + bins(1));
%     stop = Fs*(10+bins(k+1));
% 
%     AONPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
% end

time3 = Data_compilation.events.Frame(3); %odor start
time4 = Data_compilation.events.Frame(4); %odor end
time5 = Data_compilation.events.Frame(5); %post stim end

% CUMCUBEPre = AONPre;

%%

KernelType = 'Poly';
doPlot = 0;
nBoot = 3;
ANSWER_Pre = zeros(16,20,nBoot);

for odorNum = 2:16
    for boot = 1:nBoot
        odorSeq = randperm(16,odorNum);
        Results = SVM_OdorIdentifierLargeOdor_190908_AON_APC(CUMCUBEPre(:,:,odorSeq,:),KernelType,doPlot); % seems like the function looks at 4 rep so need to change that               
        ANSWER_Pre(odorNum-1,:,boot) = mean([Results.TruePos],2);  
    end
    odorNum
end
%%

figure(101);
subplot(2,1,1);hold on;
imagesc(squeeze(mean(ANSWER_Pre,3)),[0 1]);colormap('jet')
for j = 1:size(ANSWER_Pre,1)    
    figure(101);
    subplot(2,1,1);hold on;
    temp = squeeze(mean(ANSWER_Pre,3));
    idx = find(temp(j,:)>0.5);
    if ~isempty(idx)
        plot(idx(1),j,'ko')
        TP50(j,1) = idx(1);
    end
    axis([0 20 0 17])
     
end