%% Filepaths
addpath (genpath('C:\Users\Marie\Documents\Data\2021-07-23_10-27-41'))

%session 1 - conc
myKsDir = 'C:\Users\Marie\Documents\Data\2021-07-23_10-27-41'; % directory with kilosort output

%% Load data from kilosort/phy
% fct from spikes package - found in phyhelper

% sp.st are spike times in seconds (for all spikes)
% sp.clu are cluster identities (for all spikes)
% sp.cids is list of unique clusters
% sp.cgs are cluster defs (1 = MUA, 2 = good, 3 = Unsorted??) (1/cluster)
% spikes from clusters labeled "noise" have already been omitted
sp = loadKSdir(myKsDir);

%%
i = 1;
for mycluster = 1:length(sp.cids) % for each cluster 
    if sp.cgs(mycluster) == 2 | sp.cgs(mycluster) == 1  % for clusters labeled as good and mua
        %figure(i)
        %get spike times for this cluster 
        allgoodspikes = sp.st(sp.clu == sp.cids(mycluster)); % in seconds
        %t = axes;
        clusternumber = sp.cids(mycluster);
        clustertype = sp.cgs(mycluster);
        %PlotAutoCorrelogram(allgoodspikes,clusternumber,clustertype);
        PlotISI(allgoodspikes,clusternumber,clustertype)
        %acf(allgoodspikes,0.001,1);
        %i = i+1;
    end
end

function [xLin, nLin] = PlotAutoCorrelogram(st,mycluster,clustertype)
% Computes an autocorrelogram with both linear bins
% INPUT: st = spike times

%binSize = 0.0005;
binSize = 0.002;

b = 0.0001:binSize:1; % start at a tenth of a ms to avoid counting itself

[n,xLin] = histdiff(st, st, b);
nLin = n./binSize;
nLin = nLin(:); xLin = xLin(:);

refractLine = 0.002;
plotUntil = 0.1;

if clustertype == 1 % mua
    color = 'r';
elseif clustertype == 2 % good
    color = 'k';
end

b = b(:);
xx = [b(1); reshape([b(2:end-1) b(2:end-1)]', (numel(b)-2)*2,1); b(end)];
yy = reshape([nLin nLin]', numel(nLin)*2,1);
yy = yy./max(yy); % normalize 
incl = xx<plotUntil;

RefPer_duration = 0.002;
refractory_period_violation = sum(diff(st) <= RefPer_duration);
violationRate = refractory_period_violation/length(st);

plot(nexttile,xx(incl), yy(incl), 'k', 'LineWidth', 2.0);
title(['cluster', num2str(mycluster)],num2str(refractory_period_violation),'Color', color);
set(gca,'XTick', [0 0.01 0.1], 'XTickLabel', [0 10 100]);
%set(gca,'YTick', []);
box(gca,'off');
hold on;

plot(refractLine*[1 1], [0 max(yy(incl))], 'k--');
%plot([0 plotUntil], asymptVal*[1 1], 'k--');
xlim([0 plotUntil]);

if max(yy(incl))>0
    ylim([0 max(yy(incl))]);
end

end

%%
function PlotISI(st,mycluster,clustertype)

binSize = 0.002;
isi_edges = 0.0001:binSize:0.1; % start at a tenth of a ms to avoid counting itself
isi_centers = isi_edges(1:end-1)+binSize/2; % for plotting

isi = diff(st);
isih = histc(isi,isi_edges);

isih = isih./sum(isih); % normalize 
 
RefPer_duration = 0.002;
refractory_period_violation = sum(diff(st) <= RefPer_duration);
violationRate = refractory_period_violation/length(st);

refractLine = 0.002;
plotUntil = 0.1;

if clustertype == 1 % mua
    color = 'r';
elseif clustertype == 2 % good
    color = 'k';
end

bar(nexttile,isi_centers,isih(1:end-1)); % remember to ignore the last bin of histc() output
hold on;
title(['cluster', num2str(mycluster)],num2str(refractory_period_violation),'Color', color, 'FontSize', 2);
set(gca,'FontSize',20,'XLim',[0 0.25]); xlabel('ISI (s)'); ylabel('count'); grid on;
xlim([0 plotUntil]);


% plot(nexttile,xx(incl), yy(incl), 'k', 'LineWidth', 2.0);
% title(['cluster', num2str(mycluster)],num2str(refractory_period_violation),'Color', color);
% set(gca,'XTick', [0 0.01 0.1], 'XTickLabel', [0 10 100]);
% %set(gca,'YTick', []);
% box(gca,'off');
% hold on;
% 
% plot(refractLine*[1 1], [0 max(yy(incl))], 'k--');
% %plot([0 plotUntil], asymptVal*[1 1], 'k--');
% xlim([0 plotUntil]);

% if max(yy(incl))>0
%     ylim([0 max(yy(incl))]);
% end
end

