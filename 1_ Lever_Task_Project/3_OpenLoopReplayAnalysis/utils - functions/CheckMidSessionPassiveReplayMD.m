function [] = CheckMidSessionPassiveReplayMD(MySession)

% Load the relevant variables
load(MySession, 'Traces', 'TrialInfo', 'TTLs', 'ReplayTTLs', 'SingleUnits');

[OpenLoop] = ExtractReplayTrials(Traces,TrialInfo,TTLs,ReplayTTLs);

PlotUnits = 1:length(SingleUnits);

ProcessOpenLoopTrialsMD(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
    'plotfigures', 1, 'plotephys', 1, 'whichunits', PlotUnits);

end