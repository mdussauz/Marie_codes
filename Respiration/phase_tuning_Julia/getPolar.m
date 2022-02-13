function[polar] = getPolar(spike, phaseStarts) 
% Transforms a spike time into a respiration phase polar
% Input:
% - spike = a spike time
% - phaseStarts = a list of respiration phase start times
    temp = spike - phaseStarts;
    temp(temp <0) = inf;
    [~,closestIndex] = min(temp);
    polar = (spike - phaseStarts(closestIndex))/(phaseStarts(closestIndex + 1) - phaseStarts(closestIndex));
    polar = polar*2*pi;
end