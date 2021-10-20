function varargout = response_transition_points(response_key, response_time, ...
    dur_before_after, key_options)
%
% function varargout = response_transition_points(response_key, response_time, ...
%     dur_before_after, key_options)
%
% This function takes a series of continuous response keys and their time
% values and plots the distribution of durations starting with a unique
% keypress
% rem_last = last dominance duration is often incomplete so remove it
%
% Katyal, 02/2014: Wrote it
%

% First get the unique key values

key_types = unique(response_key); % Find the unique identifiers for the different states, (2 or 3 for rivalry)
if ~exist('key_options', 'var'), key_options = key_types; end
transit_points = find([true diff(response_key)~=0]); % Add a transit point at the beginning, equivalent to saying that the response at time 0 is the same as time 1

rt = response_time(transit_points);
duration_transits = [diff(rt) ...
    response_time(end)-response_time(transit_points(end))]; % Add an extra duration at the end, which is the ending of the response time-time at the last transit point
response_transits = response_key(transit_points);

num_trans = numel(transit_points);
trans_times = cell(numel(key_options));
trans_durs = cell(numel(key_options));
ndba = numel(dur_before_after);
if ndba==1, dur_before_after = [dur_before_after dur_before_after]; end

switch ndba
    case 2
        for tp = 2:num_trans
            if duration_transits(tp-1)>=dur_before_after(1) && ...
                    duration_transits(tp)>=dur_before_after(2)
                
                trans_times{response_transits(tp-1), response_transits(tp)} = ...
                    [trans_times{response_transits(tp-1), response_transits(tp)}, ...
                    rt(tp)];
                
                if nargout>1
                    trans_durs{response_transits(tp-1), response_transits(tp)} = ...
                        [trans_durs{response_transits(tp-1), response_transits(tp)}, ...
                        [duration_transits(tp-1); duration_transits(tp)]];
                end
            end
        end
        
    case 4
        for tp = 2:num_trans
            if duration_transits(tp-1)>=dur_before_after(2) && ...
                    duration_transits(tp-1)<=dur_before_after(1) && ...
                    duration_transits(tp)>=dur_before_after(3) && ...
                    duration_transits(tp)<=dur_before_after(4)
                
                trans_times{response_transits(tp-1), response_transits(tp)} = ...
                    [trans_times{response_transits(tp-1), response_transits(tp)}, ...
                    rt(tp)];
                
                if nargout>1
                    trans_durs{response_transits(tp-1), response_transits(tp)} = ...
                        [trans_durs{response_transits(tp-1), response_transits(tp)}, ...
                        duration_transits(tp)];
                end
            end
        end
end

varargout{1} = trans_times;
if nargout>1
    varargout{2} = trans_durs;
end
end