function [tt2, rr2] = interp_to_next_response(tt, rr, a, b, n, extrapval, ...
    cover_brief_mixed)
%
% [tt2, rr2] = interp_to_next_response(tt, rr, a, b, n, extrapval, ...
%    cover_brief_mixed)
%
% takes the response keys, their timings resamples and performs
% interpolation such that prior to a new point all points are filled with
% the previous value
% tt: time of response, rr: response keys, a, b: minimum and maximum values
% extrapval: initial value before response
% cover_brief_mixed: 2 element vector to interpolate short mixed percepts
% at the mid point with previous and next buton presses, [button press val,
% minimum mixed duration to allow]
% Katyal 09/2014, Wrote it
%

tt2 = linspace(a, b, n);
rr2 = zeros(size(tt2), 'single');

transit_points = find([true logical(diff(rr))]);
transit_times = [tt(transit_points) b]; % all transitions plus last time point

if ~exist('extrapval', 'var')||isempty(extrapval), extrapval = rr(1); end
if exist('cover_brief_mixed', 'var') && ~isempty(cover_brief_mixed) ...
    && cover_brief_mixed(2)>0
    cbm = true; 
else
    cbm = false;
end

rr2(tt2<transit_times(1)) = extrapval;
prevr = extrapval;

lentrm1 = length(transit_times)-1;

for tr = 1:lentrm1
    if cbm && rr(transit_points(tr))==cover_brief_mixed(1) &&...
            tr<lentrm1 && rr(transit_points(tr+1))~=cover_brief_mixed(1) &&...
            prevr~=cover_brief_mixed(1) && ...
            (transit_times(tr+1)-transit_times(tr)<=cover_brief_mixed(2))
        % mid point interpolate the brief mixed percepts
        mp = transit_times(tr) + (transit_times(tr+1)-transit_times(tr))/2;
        rr2(tt2>=transit_times(tr) & tt2<mp) = prevr;
        rr2(tt2>=mp & tt2<=transit_times(tr+1)) = rr(transit_points(tr+1));
    else
        rr2(tt2>=transit_times(tr) & tt2<=transit_times(tr+1)) = ...
            rr(transit_points(tr));
    end
    prevr = rr(transit_points(tr));
end

end