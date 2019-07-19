function bin_times = boxcar_bintimes(time_window, boxcar_size, boxcar_slide)
%bintimes = boxcar_bintimes(time_spike_vector, boxcar_size, boxcar_slide, align_time)
%
% boxcar_bintimes uses a brute force technique to find the high low and
%   mean bin times for a sliding window (boxcar) over the timewindow in
%   time_vector
%
%
% inputs
%
%   time_vector is a column vector of videosample timestamps 
%   boxcar_size is the length of the moving time window in seconds. It
%       should be divisible by the even time increments used to build
%       time_vector
%   boxcar_slide is the amount of time the boxcar should slide forward with
%       each iteration
%   align_time is the "zero time" within time_vector, and only affects the
%       output bin_center_times. If no time is input, default is the
%       earliest time in time_vector
%   time_window is the aligned window (seconds from align_time) through
%       with the sliding window will move
%
% outputs
%
%   bin_center_times is the mean between the low and high ends of the 
%       moving window as in time_vector. If there is an input align_time,
%       bin_center_times will be output as positive and negative seconds
%       from align_time.

%CHECK INPUTS 
%

    %sort time_vector
    %if ~isequal(sort(time_vector), time_vector)
    %    time_vector = sort(time_vector);
    %    warning('time_vector and spike_events sorted')
    %end

    time_vector = min(time_window):.01:max(time_window);

    if boxcar_size + boxcar_slide > (time_vector(end) - time_vector(1))
        error('boxcar properties too large for time_window')
    end

    
%CALCULATE ALL HIGH AND LOW ENDS OF BOXCAR
%
    %how many moving boxcars fit in the time window    
    bin_times = [];
    low = time_vector(1);
    high = low + boxcar_size;
    
    while high <= time_vector(end)
        
        %time_bounds = [time_bounds; [low high]];
        bin_times = [bin_times; [low mean([low high]) high]];
      
        low = low + boxcar_slide;
        high = high + boxcar_slide;
        
    end
    
%CENTER BINS AROUND ZERO
%closest to zero
    minrow = find(abs(bin_times(:,2)) == min(abs(bin_times(:,2))), 1, 'last');
    bin_times = bin_times - repmat(bin_times(minrow,2), size(bin_times));
    
end