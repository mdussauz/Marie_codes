function [Events] = ParseOpenEphysEvents(filename)

% read events and timestamps from open ephys file
[data, Timestamps, info] = load_open_ephys_data(filename);

% On how many channels were there events?
NChannels = unique(data); % data has channels IDs
% 0 is trial event / 1 is odor event 


% put On-Off timestamps for each channel in a struct
for i = 1:numel(NChannels)
   %info.eventId = 0 if off and = 1 if on
   Events.(['Channel',num2str(NChannels(i))]).On = ...
       Timestamps(intersect(find(data==NChannels(i)),find(info.eventId==1)));
   Events.(['Channel',num2str(NChannels(i))]).Off = ...
       Timestamps(intersect(find(data==NChannels(i)),find(info.eventId==0)));
end

%% for MD - what we do here: 
% C = intersect(A,B) for vectors A and B, returns the values common to
%   the two vectors with no repetitions. C will be sorted.

% Example:
% if i = 1 and eventId = 1
%
% find(data==NChannels(1) = gives all the indexes of event 0 in data
% find(info.eventId==1) = gives all the indexes of when any event was on in
% info.eventId

% So intersect here would give a vector of indices that allows to find the
% timestamps of event 0 when it went on (meaning trial ON timestamps)


