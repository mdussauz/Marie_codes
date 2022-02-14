clear;clc;close all;
dataPath = 'D:\Dropbox\MATLAB\Olfaction\MarieData\210721_dataAON_APCexperiments\NewData\LargeOdorData\';

load([dataPath 'smoothPSTH16od_AON.mat'])
temp=smoothPSTH(:,:,:,1:5);
bins = 0:0.2:4;
Fs = 1000;

AONPre = [];
for k = 1:length(bins)-1
    start = Fs*(10 + bins(1));
    stop = Fs*(10+bins(k+1));
%     start = Fs*(10+bins(k));
%     stop = start+200;
%     bsline = nanmean(temp(:,:,2000:10000,:),3);
%     intTemp = (nanmean(temp(:,:,start:stop,:),3) - bsline)./(bsline) ;
%     intTemp = (nanmean(temp(:,:,start:stop,:),3) - bsline);

    AONPre(k,:,:,:) = squeeze(nanmean(temp(:,:,start:stop,:),3));
%     AONPre(k,:,:,:) = squeeze(intTemp);
end
CUMCUBEPre = AONPre;

%%

KernelType = 'Poly';
doPlot = 0;
nBoot = 3;
ANSWER_Pre = zeros(16,20,nBoot);
ANSWER_Mus = zeros(16,20,nBoot);

for odorNum = 2:16
    for boot = 1:nBoot
        odorSeq = randperm(16,odorNum);
        Results = SVM_OdorIdentifierLargeOdor_190908_AON_APC(CUMCUBEPre(:,:,odorSeq,:),KernelType,doPlot);                
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