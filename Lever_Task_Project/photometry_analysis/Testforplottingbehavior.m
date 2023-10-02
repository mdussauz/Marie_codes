%Testfor plotting behavior
timestamps = TracesOut.Timestamps{1};
Lever = TracesOut.Lever{1};
Sniffs = TracesOut.Sniffs{1};
Licks = TracesOut.Licks{1};
Rewards = TracesOut.Licks{1};
Trial = TracesOut.Trial{1};
OdorBoxHeight = 5;
TZList =  TargetZones(TrialInfo.TargetZoneType(startTrial:endTrial),[3 1 1 3])';
TZ = TZList;
TZ = TZ';
Traces.TargetZone(perturbationTrials(thisZoneTrials(j),1)) = ...
                {TrialInfo.TargetZoneType(perturbationTrials(thisZoneTrials(j),1)) + ...
                0*foo};
%Trial(Trial~=whichOdor) = 0;
%         PlotBehaviorPerturbations(timestamps,TracesOut.Lever{1},...
%             TracesOut.Sniffs{1},TracesOut.Licks{1},TracesOut.Rewards{1},...
%             Trial,...
%             TZ, ...
%             TracesOut.HaltFlip{1},5);
        
PlotBehavior(timestamps,Lever,Sniffs,Licks,Rewards,Trial,TZ,OdorBoxHeight)