close all; clear all;

%WhichSession = 'AON_early';
%WhichSession = 'AON_late';
%WhichSession = 'APC_pre';
WhichSession = 'APC_post';

SessionPath = 'C:\Users\Marie\Documents\data\Data_Photometry_Diego';
handles.WhereSession.String = fullfile(SessionPath,WhichSession);
load(handles.WhereSession.String, 'Data_compilation')

disp('mouse x frames x odor x rep')
disp(['Dim of data is ',num2str(size(Data_compilation.Data))])

time1 = Data_compilation.events.Frame(1);
time2 = Data_compilation.events.Frame(2);
time3 = Data_compilation.events.Frame(3);
time4 = Data_compilation.events.Frame(4);
time5 = Data_compilation.events.Frame(5);

%% Plot raw data 
for mouse = 1:size(Data_compilation.Data,1)
    figure(mouse);
    sgtitle(['mouse #', num2str(mouse)])
    for odor = 1:16
        subplot(4,4,odor)
        for rep = 1:size(Data_compilation.Data,4)
            plot(Data_compilation.Data(mouse,:,odor,rep)); hold on;
        end
        xline([time1, time2, time3, time4, time5]);
    end
end

%% Plot mean and variance across repeats 
for mouse = 1:size(Data_compilation.Data,1)
    figure(mouse+size(Data_compilation.Data,1));
    sgtitle(['mouse #', num2str(mouse)])
    for odor = 1:16
        x = 1:size(squeeze(Data_compilation.Data(mouse,:,odor,:)),1);
        y = squeeze(Data_compilation.Data(mouse,:,odor,:));
        if sum(y) ~= 0
            odor_mean = mean(y,2);
            odor_std = std(y,0,2);
            subplot(4,4,odor)
            %shadedErrorBar(x,odor_mean,odor_std,'lineProps','g');
            shadedErrorBar(x,y',{@mean, @std},'lineProps','g');
            xline([time1, time2, time3, time4, time5]);
        else
            disp(['no activity for mouse', num2str(mouse)])
        end
    end
end