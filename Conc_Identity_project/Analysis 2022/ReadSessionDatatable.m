function SessionInfo = ReadSessionDatatable(FilePath, FileName)
% -- written by MD
% -- read excel table containing session info 
% return a structure with the following fields: 
% RecFile: name of recording file (date + time)
% StimFile: name of stimulus file containing info about the trials
% SessionType: conc or id, ie concentration or identity experiment
% MouseName: name of subject
% BrainRegion: recorded in AON or APC

if nargin < 1
    if strcmp(computer,'PCWIN64')
        FilePath = 'C:\Users\Marie\Documents\Code\Marie_codes\Conc_Identity_project\Analysis 2022';
    else
        FilePath = 'opt/Marie_codes/Conc_Identity_project/Analysis 2022';
    end
  
  FileName = 'Final_Table_Conc_Id_Analysis.xlsx'; 
end

T = readtable(fullfile(FilePath,FileName));
T = table2struct (T,"ToScalar",true);

SessionInfo = T;

end
