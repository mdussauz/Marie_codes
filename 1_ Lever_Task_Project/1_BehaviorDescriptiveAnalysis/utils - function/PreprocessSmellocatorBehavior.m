function [TrialInfo, Traces, TargetZones] = PreprocessSmellocatorBehavior(MyFilePath)
% Simplified version of PreprocessSmellocatorData to get cleaned up behavior info 
% It parses behavioral data from the mouse lever task
% into trials, with relevant continuous (lever, motor, respiration, lickpiezo)
% and event data (licks, target zone flags, odor ON-OFF, etc) for each trial

%% globals
global SampleRate;
SampleRate = 500; % Samples/second
global startoffset;
startoffset = 1; % in seconds
global savereplayfigs;
savereplayfigs = 0;
global errorflags; % [digital-analog sample drops, timestamp drops, RE voltage drift, motor slips]
errorflags = [0 0 0 0];
global TargetZones; 

%% core data extraction (and settings)
[MyData, MySettings, DataTags] = ReadSessionData(MyFilePath);
FileLocations.Behavior = MyFilePath;
[FilePaths, MyFileName] = fileparts(MyFilePath);
disp(MyFileName);

%% Parse into trials
[Trials,InitiationsFixed] = CorrectMatlabSampleDrops(MyData, MySettings, DataTags);
[MyData, DataTags] = OdorLocationSanityCheck(MyData, DataTags); % check that manifold actually moved as expected
[Traces, TrialInfo] = ParseBehavior2Trials(MyData, MySettings, DataTags, Trials);

% sanity check - did some guess work in CorrectMatlabSampleDrops to compute
% odor start - check if it made sense
if ~isempty(InitiationsFixed)
    if any(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01)
        weirdo = find(abs(diff(TrialInfo.OdorStart(InitiationsFixed,:),1,2))>=0.01);
        if any(TrialInfo.OdorStart(InitiationsFixed(weirdo),2)>-1)
            disp('something funky with computing odorstart from Lever trace');
            %keyboard;
            TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
        else
            % Initiation hold was larger than a second - that's couldn't
            % compute it accurately from trial traces
            TrialInfo.OdorStart(InitiationsFixed(weirdo),1) = TrialInfo.OdorStart(InitiationsFixed(weirdo),2);
        end
    end
end

%% Check if passive tuning was done
MyTuningTrials = []; TuningTrialSequence = []; PassiveReplayTraces = []; TuningParams = [];
[TuningFile] = WhereTuningFile(FilePaths,MyFileName);
if ~isempty(TuningFile)
    [MyTuningTrials, TuningTrialSequence, PassiveReplayTraces, TuningParams] = ParseTuningSession(TuningFile);
    disp(['Found Tuning File: ',TuningFile]);
    FileLocations.Tuning = TuningFile;
end

%% Sniffs - ERROR HERE TO BE SOLVED
% [SniffTS] = ReadThermistorData(MyFilePath); % in behavior timestamps
% [SniffTS_passive] = ReadThermistorData(TuningFile); % in behavior timestamps


%% clearing up disk space
clearvars -except TrialInfo Traces TargetZones ...  
    PassiveReplayTraces  startoffset SampleRate ...
   % SniffTS SniffTS_passive
end