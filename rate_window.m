function [rate_mtx] = rate_window(eptrials, lost_trials, clusters, wndw, bins, boxcar_shift) %flag_events
%finds the rate of every cell during each bin on every trial in the window
%around flag
%
%window is a 1x2 row matrix with the beinging and end of the time window
%expressed in seconds relative to the flag events. e.g. [-1 1] is one
%second before to one second after
%
%bin size is the size of the sliding window used to calculate firing rates
%across the window


%correct trials
corrects = unique(eptrials(eptrials(:,9)==1, 6));
corrects = setdiff(corrects, lost_trials);

%errors = unique(eptrials(eptrials(:,9)==0, 6));

%flag times
flag_events = nan(size(corrects));
for trl_num = 1:length(corrects)

    %event time
    %et = min(eptrials(eptrials(:,6)==trl_num & ismember(eptrials(:,13), [1 2]), 1)); %lighton
    et = min(eptrials(eptrials(:,6)==trl_num, 1)); %trial start
    
    if ~isempty(et)
        %for now, light on
        flag_events(trl_num) = et;
    end
end

%calculate rates over moving window
sw_bins = floor(sum(abs(wndw))/boxcar_shift - (bins/boxcar_shift - 1)); %how many sliding windows we can fit
rate_mtx = nan(length(clusters), sw_bins, length(corrects));



for trl_num = 1:length(corrects)
    
    %all low ends of the sliding windows
    common_low = (0:boxcar_shift:((sw_bins-1)*boxcar_shift));
    bins_low = common_low + repmat(flag_events(trl_num) + wndw(1), 1, sw_bins);

    %calculate spike counts for every cluster at each sliding window point
    %over the window
    for bin_n = 1:length(bins_low)
        low_idx = eptrials(:,1) >= bins_low(bin_n);
        high_idx = eptrials(:,1) <= bins_low(bin_n)+bins;
        
        rate_mtx(:,bin_n, trl_num) = histc(eptrials(low_idx & high_idx, 4), clusters);

    end
end







end