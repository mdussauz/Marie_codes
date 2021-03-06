%% OVERVIEW
%{

1 - load some data from one tetrode 
2 - remove very low  and very high frequency data, leaving only the freq.
band where spikes live
3 - apply a threshold to find the times of these spikes
4- for each of these spikes, remember the max. amplitude, and the waveform
5 - look at the clustering of amplitudes across tetrode channels

%}


%% 1- load raw data

set(0,'DefaultFigureWindowStyle','docked'); % fix matlab's figure positioning bug

% raw data available on
% https://drive.google.com/drive/folders/1CwFcErgp3F3D6I2TB_hTtW1JAQB21TAC?usp=sharing
%
%datapath='C:\Users\mdussauz\Desktop\Spike_Sorting\TENSS\spike_sorting_tutorial_Nacho\tenss_example_data'; % <- edit this
datapath='/mnt/data/N5/2019-09-24_15-13-13'

data_raw=[];
channelid = [1, 29:31];
for ch=1:4 % grab 4 channels of raw data from one tetrode
    fname = sprintf('100_CH%d.continuous',channelid(ch));
    fprintf('loading ch %d/4 [%d]\n',ch,channelid(ch));
    [data, timestamps, info] = load_open_ephys_data_faster(fullfile(datapath,fname));
   % data=data(1:20000); % cut down data
    data_raw(:,ch) = data;
end
disp('done');

size(data_raw);

data_raw = data_raw.*info.header.bitVolts;
fs = info.header.sampleRate;

%% plot some raw data

plotlim=20000;

% 1 electrode (only run the next 2 lines of code for now)
figure(1);
plot(data_raw(1:plotlim,:));

% tetrodes (now only run the next two lines of code)
figure(2);

plot(data_raw(1:plotlim,:)+100*repmat([1:4]',1,plotlim)');



%% filter exercise 1: change freq band


[b1,a1] = butter(1, [300 6000]/(fs/2),'bandpass'); % filter 1 (normalize bp freq. to nyquist freq.)
data_bp1=filter(b1,a1,data_raw(1:plotlim,:)); % apply filter 1 
figure(3); clf;
plot(data_bp1(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','k');


[b2,a2] = butter(1, [300 500]/(fs/2),'bandpass'); % filter 2
data_bp2=filter(b2,a2,data_raw(1:plotlim,:)); % apply filter 2 
hold on
plot(data_bp2(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','b');

%% filter exercise 2: change filter order

[b1,a1] = butter(1, [100 6000]/(fs/2),'bandpass'); % choose filter (normalize bp freq. to nyquist freq.)
[b2,a2] = butter(7, [100 6000]/(fs/2),'bandpass');

data_bp1=filter(b1,a1,data_raw(1:plotlim,:)); %apply filter 1
data_bp2=filter(b2,a2,data_raw(1:plotlim,:)); %apply filter 2

figure(4);
p1=plot(data_bp1(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','k');hold on
p2=plot(data_bp2(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','r');
legend([p1(1),p2(1)],'first order','fifth order')


% Q: why not always use the highest order?
%% filter exercise 3: causal vs. acuasal


[b1,a1] = butter(4, [300 6000]/(fs/2)); % choose filter (normalize bp freq. to nyquist freq.)

data_bp_filt=filter(b1,a1,data_raw(1:plotlim,:)); %apply filter 1 in one direction
data_bp_filtfilt=filtfilt(b1,a1,data_raw(1:plotlim,:)); %apply filter 1 in both directions

figure(5);
p1=plot(data_bp_filt(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','k-');hold on
p2=plot(data_bp_filtfilt(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','r-');hold on
p3=plot(data_raw(1:plotlim,:)+100*repmat([1:4]',1,plotlim)','g-');
legend([p1(1),p2(1),p3(1)],'filter','filtfilt','raw')
%% 2 -pick your filter and apply for spike sorting

% duration_to_anlyse = 5*60*fs;
duration_to_anlyse = size(data_raw, 1);

[b,a] = butter(1,[300 6000]/(fs/2),'bandpass'); 
data_bp = filtfilt(b,a,data_raw(1:duration_to_anlyse,:));

%% 3 - thresholding
% find treshold crossings
threshold=-15;
crossed= min(data_bp,[],2) < threshold; % trigger if _any_ channel crosses in neg. direction

spike_onsets=find(diff(crossed)==1);

% figure;plot(data_bp(1:10000,1),'k');hold on;plot(10*crossed(1:10000),'r')

length_sec=duration_to_anlyse/fs;
fprintf('got %d candidate events in %dmin of data, ~%.2f Hz\n',numel(spike_onsets),round(length_sec/60),numel(spike_onsets)/length_sec);

figure(6);
plotlim = 200000;
for ch=1:4;
    for i=1:numel(spike_onsets)
        if(spike_onsets(i)<plotlim)
            plot([1 1].*spike_onsets(i),100*ch+[-1 1].*threshold*2,'k-');hold on
        end
    end
    plot(100*ch+data_bp(1:plotlim,ch));
end
%  change threshold

%% 4 - extract spike waveforms and make some features

spike_window=[1:32]-5; % the -5 is there to grab some pre-treshold crossing samples

spikes=[];
spikes.waveforms=zeros(numel(spike_onsets),4*numel(spike_window)); % pre-allocate memory
spikes.peakamps=zeros(numel(spike_onsets),4);
spikes.times = spike_onsets/(fs/1000);

for i=1:numel(spike_onsets)
    this_spike=(data_bp(spike_onsets(i)+spike_window,:));
    
    spikes.waveforms(i,:)= this_spike(:);% grab entire waveform
    spikes.peakamps(i,:)=min(this_spike); % grab 4 peak amplitudes
end


%% 5 - plot peak to peak amplitudes
figure(7);
plot(spikes.peakamps(:,1),spikes.peakamps(:,2),'.'); % try a few different combinations, like 1 vs. 2, 2 vs. 3 etc.
daspect([1 1 1]);

%% initialize all cluster assignments to 0
spikes.cluster=zeros(numel(spike_onsets),1);

%% manual spike sorter
% % This code now makes use of some fun tricks - dont feel like you need to
% % understand all of this.
% %
% %
% % cluster 0 shall be the noise cluster (we dont plot this one)
% 
% % Commands:
% % up/down : cycle through projections
% % 0-9: Select cluster to examine and change
% % + : Draw contour to add to currently selected cluster
% % * : Draw polygon and add all _non-selected_ spikes to noise cluster (0)
% % t : toggle display of currently selected cluster
% 
% 
% figure(1);
% run =1;
% 
% projections=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4]; % possible feature projections
% use_projection = 1;
% 
% cluster_selected = 2;  
% spike_selected = 1;
% show_cluster = true(1,10);
% cluster_colors = [[0,0,0]; hsv(9)];
% 
% while run
%     
%     dat_x = spikes.peakamps(:,projections(use_projection,1));
%     dat_y = spikes.peakamps(:,projections(use_projection,2));
%     
%     clf;
%     subplot(2,3,1); hold on; % plot median waveform
%     plot(quantile(spikes.waveforms(spikes.cluster==cluster_selected,:),.5),'k', 'Linewidth',2);
%     plot(quantile(spikes.waveforms(spikes.cluster==cluster_selected,:),.8),'k--');
%     plot(quantile(spikes.waveforms(spikes.cluster==cluster_selected,:),.2),'k--');
%     
%     
%     plot(spikes.waveforms(spike_selected,:),'r', 'Linewidth',2); % also plot currently selected spike waveform
%     
%     title('waveforms from cluster');
%     
%     subplot(2,3,4); hold on; %  plot isi distribution
%     isi = diff(spikes.times(spikes.cluster==cluster_selected));
%     bins = linspace(0,100,50);
%     h= hist(isi,bins); h(end)=0;
%     stairs(bins,h, 'k-', 'linewidth',2);
%     title('ISI histogram'); xlabel('isi(ms)');
%     
%     ax=subplot(2,3,[2 3 5 6]); hold on; % plot main feature display
%     
%     %b_show = ones(size(spikes.cluster));
%     %ii = spikes.cluster > 0; %  plot noise cluster
%     ii = true(size(spikes.cluster));
%     for i = 1:numel(show_cluster) % remove hidden units
%         if ~show_cluster(i)
%            ii(spikes.cluster == i - 1)  =  false;
%         end
%     end
%     scatter(dat_x(ii), dat_y(ii), (0.5+(spikes.cluster(ii)==cluster_selected))*20, cluster_colors(spikes.cluster(ii)+1, :), 'filled');
%     plot(dat_x(spike_selected),dat_y(spike_selected),'ro','markerSize',10);
%     xl = xlim(gca);
%     dot_x = linspace(xl(1), xl(end),10);
%     dot_y = zeros(size(dot_x));
%     scatter(dot_x(show_cluster), dot_y(show_cluster), 100 * ones(size(dot_x(show_cluster))), cluster_colors(show_cluster, :), 'filled','d');
%     title(sprintf('current cluster %d, projection %d, %d spikes in cluster', cluster_selected, use_projection, sum(spikes.cluster==cluster_selected)));
%     xlabel('0-9:select cluster, +: add to cluster, -: remove,  up/down: change projection');
%     grid on;
%     
%     [x,y,b]=ginput(1);
%     
%     if b>47 & b <58 % number keys, cluster select
%         cluster_selected=b-48;
%     end
%     
%     % Toggle showing a cluster
%     if b == 116
%         if show_cluster(cluster_selected + 1) 
%              show_cluster(cluster_selected + 1)  = false;
%         else
%              show_cluster(cluster_selected + 1)  = true;
%         end
%    
%     end
%     if b==30; use_projection=mod(use_projection,6)+1; end; % up/down: cycle trough projections
%     if b==31; use_projection=mod(use_projection-2,6)+1; end; % up/down: cycle trough projections
%     if b==27; disp('exited'); run=0; end; % esc: exit
%     
%     if b==43 | b==42 | b==95; % +, add to cluster
%         t= imfreehand(ax,'Closed' ,1);
%         t.setClosed(1);
%         r=t.getPosition;
%         px=r(:,1);py=r(:,2);
%         in = inpolygon(dat_x,dat_y,px,py);
%         if b==43 % +, add
%             spikes.cluster(in)=cluster_selected;
%         elseif b==95 % -: subtract from cluster, to noise
%             spikes.cluster(in)=0;
%         else % *. intersect cluster (move all non selected to null cluster)
%             spikes.cluster(~in & spikes.cluster==cluster_selected)=1;
%         end;
%     end;
%     
%     if b==1 % left click - select individual waveform to plot
%         [~,spike_selected]=min((dat_x-x).^2 +(dat_y-y).^2);
%     end;
%     
% end;
