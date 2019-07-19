function [counts, absolute_times] = spike_count_boxcar(time_spike_vector, clusters, bin_times, align)
%[counts, bin_center_times] = spike_count_boxcar(time_spike_vector, clusters, boxcar_size, boxcar_slide, align_time)
%
% spike_count_boxcar counts the number of spikes from each of the cells in
% 'clusters' over times in 'time_vector.' It calculates spike counts
% using a sliding boxcar method whereby a sliding window of length 
% 'boxcar_size' second iterates across the time window, progressing 
% 'boxcar_slide' seconds each iteration. At each iteration, the number of
% items in spike events, indexed to include on the relevant times in
% time_vector, that also appear in clusters are counted. The center of the
% indexed time_vector times is recorded as the time for that iteration.
% If there is an input 'align_time', the center times are recorded as a
% difference from 'align_time.'
%
% inputs
%
%   time_spike_vector is a 2 by n matrix with timestamps in column 1 and 
%   clusters IDs indicating spike events in column 2. Column 2 rows without 
%       spike events should be 1s.
%   clusters is unique(spike_events) or a subset of these values indicating
%       which spike events should be counted. MAKE THIS 1 FOR MOVING
%       AVERAGE (e.g., for inst veloc or accel)
%   boxcar_size is the length of the moving time window in seconds. It
%       should be divisible by the even time increments used to build
%       time_vector
%   boxcar_slide is the amount of time the boxcar should slide forward with
%       each iteration
%   align_time is the "zero time" within time_vector, and only affects the
%       output bin_center_times. If no time is input, default is the
%       earliest time in time_vector
%
% outputs
%
%   count is a matrix of spike counts for every cell in 'clusters' over the
%       time window. It is not as long as time_vector due to the size of
%       the boxcar. ALTERNATIVELY, MOVING AVERAGE
%   bin_center_times is the mean between the low and high ends of the 
%       moving window as in time_vector. If there is an input align_time,
%       bin_center_times will be output as positive and negative seconds
%       from align_time.
%    
%


%CHECK INPUTS 
%
    %sort time_vector
    if ~isequal(sort(time_spike_vector(:,1)), time_spike_vector(:,1))
        time_spike_vector = sortrows(time_spike_vector,1);
        warning('time_vector and spike_events sorted')
    end
    
%CONSTRAIN TIME
%
time_spike_vector = time_spike_vector(time_spike_vector(:,1)>=(align+bin_times(1, 1)) & time_spike_vector(:,1)<=(align+bin_times(end, 3)), :);
   
    
%CHECK CLUSTERS
if isempty(ismember(unique(time_spike_vector(time_spike_vector(:,2)~=1,2)), clusters))
    error('spike_events and clusters must share ~nan items')
    
elseif isequal(clusters,1)
    %CALCULATE MOVING AVERAGE OF COLUMN VALUE
    %
    counts = nan(size(clusters,1), size(bin_times, 1));
    absolute_times = nan(size(bin_times, 1), 2);
    for i = 1:size(bin_times, 1)
        
        low = align + bin_times(i, 1);
        high = align + bin_times(i, 3);
        
        absolute_times(i, :) = [low high];
        
        %find time bounds
        low_idx = time_spike_vector(:,1) >= low;
        high_idx = time_spike_vector(:,1) <= high;
        
        %count
        counts(:,i) = nanmean(time_spike_vector(low_idx & high_idx, 2));
    end 
    
    
else
    %COUNT SPIKE EVENTS WITHIN BOUNDS OF bin_times
    %
    counts = nan(size(clusters,1), size(bin_times, 1));
    absolute_times = nan(size(bin_times, 1), 2);
    for i = 1:size(bin_times, 1)
        
        low = align + bin_times(i, 1);
        high = align + bin_times(i, 3);
        
        absolute_times(i, :) = [low high];
        
        %find time bounds
        low_idx = time_spike_vector(:,1) >= low;
        high_idx = time_spike_vector(:,1) <= high;
        
        %count
        counts(:,i) = histc(time_spike_vector(low_idx & high_idx, 2), clusters);
    end 
end

 

end