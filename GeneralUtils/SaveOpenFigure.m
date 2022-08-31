%% Save Open Figures

% Enter information
AnimalName = 'O7';
SessionDate = '20220630';
SessionType = 'OL_Tuning'; %OL %OL_Tuning

% Create filename 
figureFolder = '/mnt/data/Lever_Task_Figures';

%savepath = '/mnt/data/Lever_Task_Figures/O2_20211011.pdf';

savepath = fullfile(figureFolder,AnimalName,[AnimalName,'_',SessionDate,'_',SessionType,'.pdf']);

% Check if file already exists
if exist(savepath)
    reply = input('These figures have already been saved. \nDo you want to overwrite? Y/N [Y]: ','s');
    if strcmp(reply,'N')
        return;
    end
end

% Get handles to all open figures
figHandles = findall(0,'Type','figure'); 

% Save first figure
exportgraphics(figHandles(1),savepath)

% Loop through figures 2:end
for i = 2:numel(figHandles)
    exportgraphics(figHandles(i),savepath,'Append',true)
end

close all; clear all;
