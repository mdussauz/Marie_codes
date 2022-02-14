FileName = 'J1_20200722_o0.mat';

%% load the file
Temp = load(FileName,'session_data');
MyTraces = Temp.session_data.trace;
MyParams = Temp.session_data.params;
MyTrials = Temp.session_data.TrialSequence;
[nrows, ~] = size(MyTraces);


% [Lever RotaryEncoder TrialON InTargetZone InRewardZone Rewards
% Licks HomeSensor Camera1 Camera2]
MyData = MyTraces(:,[1 3 7:11]);

% add 3 cols in the beginning [timestamps TZoneHighLim TZoneLowLim ...
% add 2 cols in the end [FZoneHighLim FZoneLowLim]
MyData = horzcat(Temp.session_data.timestamps, zeros(nrows,2), MyData, zeros(nrows,2));

% append motor location
MyData(:,13) = MyTraces(:,4);

if find(ismember(Temp.session_data.trace_legend,'homesensor'))
    whichcol = find(ismember(Temp.session_data.trace_legend,'homesensor'));
    MyData(:,14) = MyTraces(:,whichcol);
end

if find(ismember(Temp.session_data.trace_legend,'respiration'))
    whichcol = find(ismember(Temp.session_data.trace_legend,'respiration'));
    MyData(:,15) = MyTraces(:,whichcol);
end

if find(ismember(Temp.session_data.trace_legend,'camerasync'))
    whichcol = find(ismember(Temp.session_data.trace_legend,'camerasync'));
    MyData(:,16) = MyTraces(:,whichcol);
end

if size(MyTraces,2)>whichcol % if there was a second camera
    whichcol = whichcol + 1;
    MyData(:,17) = MyTraces(:,whichcol);
end

clear Temp
