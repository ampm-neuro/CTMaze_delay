function [out] = common_bin_time(wndw, bins, boxcar_shift)
%high and low averaging window of each bin

    %how many sliding windows we can fit?
    sw_bins = floor(sum(abs(wndw))/boxcar_shift - (bins/boxcar_shift - 1)); 
    
    %all low ends of the sliding windows
    bins_low = (0:boxcar_shift:((sw_bins-1)*boxcar_shift)) + repmat(wndw(1), size((0:boxcar_shift:((sw_bins-1)*boxcar_shift))));

    
    out = nan(length(bins_low), 2);

    %calculate spike counts for every cluster at each sliding window point
    %over the window
    for bin_n = 1:length(bins_low)
        
        out(bin_n, :) = [bins_low(bin_n) bins_low(bin_n)+bins];

       
    end

%plot
    low = find(round(abs(out(:,2)),3) == round(min(abs(out(:,2))),3), 1, 'first')
    mid = mean(find(round(abs(mean(out,2)),3) == round(min(abs(mean(out,2))),3)))
    high = find(round(abs(out(:,1)),3) == round(min(abs(out(:,1))),3), 1, 'last')

    hold on; plot([low low], [0 .6], 'r-')
    hold on; plot([mid mid], [0 .6], 'r-')
    hold on; plot([high high], [0 .6], 'r-')
    
    
end