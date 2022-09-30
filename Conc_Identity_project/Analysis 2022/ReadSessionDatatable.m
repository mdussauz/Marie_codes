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
  FilePath = 'C:\Users\Marie\Desktop';
  FileName = 'Final_Table_Conc_Id_Analysis.xlsx'; 
end

T = readtable(fullfile(FilePath,FileName));
T = table2struct (T,"ToScalar",true);

SessionInfo = T;

end