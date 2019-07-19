function [outlying_trials outlier_index] = stem_outliers(eptrials, lost_trials)
%identifies trials where the rat spent an unusual amount of time in the stem 

    %preallocate
    stem_time = nan(max(eptrials(:,6)),1);
    
    %iterate for stem time
    for trial = setdiff(1:max(eptrials(:,6)), lost_trials)
        
        stem_ent = min(eptrials(eptrials(:,6) == trial & ismember(eptrials(:,11), 1:4), 1));
        stem_ext = max(eptrials(eptrials(:,6) == trial & ismember(eptrials(:,11), 1:4), 1));
        
        stem_time(trial) = stem_ext - stem_ent;
     
    end

    %stem_time
    median_time = nanmedian(stem_time);
    std_time = nanstd(stem_time);
    differences = abs(stem_time - repmat(median_time, size(stem_time)));
    
    %outlier index (at least 1s on stem and more than 2sds from mean)
    outlier_index = differences > std_time*2 & stem_time > 1;
    outlier_index(lost_trials) = 1;
    
    %outlying trials
    trials = 1:max(eptrials(:,6));
    outlying_trials = trials(outlier_index);

end