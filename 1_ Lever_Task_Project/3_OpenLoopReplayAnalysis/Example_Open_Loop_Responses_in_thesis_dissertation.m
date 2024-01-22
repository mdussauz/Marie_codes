%Example_Open_Loop_Responses_in_thesis_dissertation

MySession = '/Users/mariedussauze/Desktop/Analysis/data/Smellocator/Processed/Behavior/O3/O3_20211005_r0_processed.mat';

%% Select Unit to plot

unitsperfig = 4;

%Modulated units 
%unit 50 = tet9.1 / unit 15 = tet 5.1 / unit 4 = tet 1 / unit 33 = tet 7 > unit 11 to 14 on figure

%Unmodulated units
% unit 48 = tet 9 / unit 14 = tet 5.1 / unit 1 = tet 1 / unit 34 = tet 7 > unit 7 to 10 on figure
MyUnits = [48 14 1 34 50 15 4 33]; 

%% get the processed data loaded
load(MySession, 'Traces', 'PassiveReplayTraces', 'TrialInfo', ...
                'SingleUnits', 'TTLs', 'ReplayTTLs', 'TuningTTLs', ...
                'SampleRate', 'startoffset', 'TargetZones', 'errorflags');
            
OpenLoop = ExtractReplayTrials(Traces, TrialInfo, TTLs, ReplayTTLs);

ProcessOpenLoopTrials(OpenLoop, TrialInfo, SingleUnits, TTLs, ...
        'PSTHsmooth', 100, ...
        'plotfigures', 0, 'plotephys', 1, 'UnitsPerFig', unitsperfig, ...
        'whichunits', MyUnits);

