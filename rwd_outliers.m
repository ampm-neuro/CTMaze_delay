function [outlying_trials outlier_index] = rwd_outliers(eptrials, lost_trials)
%identifies trials where the rwd flag is very far from the average flagged
%loc
    
    %reward average pos
    L_rwd_pos = eptrials(~ismember(eptrials(:,6), lost_trials) & eptrials(:,5)==1, [2 3 6]); %get trial numbers
    R_rwd_pos = eptrials(~ismember(eptrials(:,6), lost_trials) & eptrials(:,5)==2, [2 3 6]);
    L_rwd_loc = nanmean(L_rwd_pos(:,1:2));
    R_rwd_loc = nanmean(R_rwd_pos(:,1:2));
    trial_nums = [L_rwd_pos(:,3); R_rwd_pos(:,3)];
    
    %preallocate
    L_dists = nan(size(L_rwd_pos,1),1);
    R_dists = nan(size(R_rwd_pos,1),1);
    
    %calculate distances
    for L = 1:length(L_dists)

        L_dists(L) = pdist([L_rwd_pos(L,1:2); L_rwd_loc]);
        
    end
    
    for R = 1:length(R_dists)
     
        R_dists(R) = pdist([R_rwd_pos(R,1:2); R_rwd_loc]);
        
    end
    
    median_L_dist = median(L_dists);
    median_R_dist = median(R_dists);
    std_L_dist = std(L_dists);
    std_R_dist = std(R_dists);
    L_differences = abs(L_dists - repmat(median_L_dist, size(L_dists)));
    R_differences = abs(R_dists - repmat(median_R_dist, size(R_dists)));
    
    %outlier indices (dist of at least 25 and more than 2sds from mean)
    L_outlier_index = L_differences > std_L_dist*2 & L_dists > 15;
    R_outlier_index = R_differences > std_R_dist*2 & R_dists > 15;

    outlying_trials = trial_nums(logical([L_outlier_index;R_outlier_index]));


end