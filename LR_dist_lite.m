function [dist_mtx, stem_ent,stem_ext, rwd_time]  = LR_dist_lite(eptrials, clusts, incl_trials, tw)
%function [] = LR_dist_lite()
%
%LR_dist finds the distance between the instantaneous population rep and
%the session representation of that session from each L and R trials and then
%outputs the difference between them
%
%tw should be a multiple of 10ms
%
%figure;


%test cell and trial proportions
num_Ltrials = sum(ismember(unique(eptrials(eptrials(:,7)==1, 6)), incl_trials));
num_Rtrials = sum(ismember(unique(eptrials(eptrials(:,7)==2, 6)), incl_trials));

if num_Ltrials-1 <= length(clusts)
    clusts = clusts(1:num_Ltrials-1)
end
if num_Rtrials-1 <= length(clusts)
    clusts = clusts(1:num_Rtrials-1)
end

%arbitrary for making dist_mtx big enough to hold any trial
arb = 25000;

%all trials in session
all_trls = unique(eptrials(:,6));

%preallocate trl_rate_mtx
trl_rate_mtx = nan(length(all_trls), length(clusts));
trial_type = nan(size(all_trls));
trial_length = nan(size(all_trls));



%mean rate for each cell on each trial type
    for trl = all_trls'
        
        %skip bad trials and errors
        if ~ismember(trl, incl_trials)
            continue
        end

        %trial_type
        trial_type(trl) = mode(eptrials(eptrials(:,6) == trl, 7));

        %4 seconds after reward
        %trl_end = max(eptrials(eptrials(:,6)==trl, 1));
        trl_end = min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,5), [1 2 3 4]), 1)) + 4;
        
        %number of vid sample in trial
        trial_length(trl) = sum(eptrials(eptrials(:,6)==trl & eptrials(:,4)==1, 1));
        trl_time = trl_end - min(eptrials(eptrials(:,6)==trl, 1));
        
        
        for c = 1:length(clusts)
            %fill rates
            trl_rate_mtx(trl, c) = length(eptrials(eptrials(:,1) <= trl_end & eptrials(:,6) == trl & eptrials(:,4) == clusts(c), 1))/trl_time*tw;
        end
    end    
    
%preallocate matrix (TW ms sliding windows in longest trial X trials)
vid_samples = max(trial_length);
bins_in_tw = tw/.01;
wdw_bins = vid_samples - bins_in_tw + 1;
dist_mtx = nan(length(all_trls), arb*2);
time_set_align = nan(length(all_trls), int32(wdw_bins));

light_on_time = nan(length(all_trls), 1);
stem_ent = nan(length(all_trls), 1);
stem_ext = nan(length(all_trls), 1);
rwd_time = nan(length(all_trls), 1);


%for all trials
    for trl = all_trls'
        
        %trl
        
        %skip bad trials and errors
        if ~ismember(trl, incl_trials)
            continue
        end
        
        %1 second after reward
        trl_end = min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,5), [1 2 3 4]), 1)) + 1;
        
        %mean for L and R excluding current trial
        L_rep = trl_rate_mtx(ismember(all_trls, incl_trials) & trial_type==1 & all_trls~=(trl),:);
        R_rep = trl_rate_mtx(ismember(all_trls, incl_trials) & trial_type==2 & all_trls~=(trl),:);
        
        %align
        light_on_time(trl) = min(eptrials(eptrials(:,6)==trl & ismember(eptrials(:,13),[1 2]),1));
        
        %min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,11), 1:4), 1)) - light_on_time(trl)
        
        stem_ent(trl) = min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,11), 1:4), 1)) - light_on_time(trl);
        stem_ext(trl) = max(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,11), 1:4), 1)) - light_on_time(trl);
        choice_exit = max(eptrials(eptrials(:,6) == trl & eptrials(:,10) == 3 & eptrials(:,1) < min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,5), [1 2 3 4]), 1)), 1));
        rwd_time(trl) = min(eptrials(eptrials(:,6) == trl & ismember(eptrials(:,5), [1 2 3 4]), 1)) - light_on_time(trl);
        
        %all available time points given window size
        time_set = eptrials(eptrials(:,1) <= trl_end & eptrials(:,6)==trl & eptrials(:,4)==1,1);
        
        %set light on time to zero
        time_set_align_hold = time_set - repmat(light_on_time(trl), size(time_set));        
            %if even number of bins in time window
            if rem(bins_in_tw, 2) == 1
                time_set_align_hold = time_set_align_hold(floor(bins_in_tw/2) : length(time_set) - floor(bins_in_tw/2));
                time_set_align(trl, floor(bins_in_tw/2) : length(time_set) - floor(bins_in_tw/2)) = time_set_align_hold;
            else

                time_set_align_hold = time_set_align_hold(1 + floor(bins_in_tw/2) : length(time_set) - floor(bins_in_tw/2));
 
                time_set_align(trl, 1 + floor(bins_in_tw/2) : length(time_set) - floor(bins_in_tw/2)) = time_set_align_hold;
            end

        %find dist between inst and L and R using sliding TW ms time window 
        %along entire trial (probably fill a cell)
        %
        %for time window
        inst_rep = nan(length(time_set) - bins_in_tw + 1, length(clusts));
        
        for i = 1:(length(time_set) - bins_in_tw + 1)
            time_point_low = time_set(i);
            time_point_hi = time_set(i+bins_in_tw-1);
            
                %find instant representation
                inst_rep(i, :) = histc(eptrials(eptrials(:,1)>=time_point_low & eptrials(:,1)<=time_point_hi, 4), clusts)';
        end

        %mahal(L_rep, R_rep')/sqrt(length(clusts))
        L_dist = sqrt(mahal(inst_rep, L_rep))/sqrt(length(clusts));
        R_dist = sqrt(mahal(inst_rep, R_rep))/sqrt(length(clusts));
        
        %standardize
        %
        L_dist = L_dist - repmat(mean(L_dist), size(L_dist));
        L_dist = L_dist./std(L_dist);
        R_dist = R_dist - repmat(mean(R_dist), size(R_dist));
        R_dist = R_dist./std(R_dist);
        %}
        
         %plot(time_set_align, smooth(L_dist, 50))
        %hold on; plot(time_set_align, smooth(R_dist, 50), 'r')

        %
        %figure; 
        
        %if left trial
        if trial_type(trl) == 1
            
            %higher values when left is smaller (inst rep is similar to left)
            dist_mtx(trl, 1:length(L_dist)) = R_dist - L_dist;
            %dist_mtx(trl, 1:length(L_dist)) = dist_mtx(trl, 1:length(L_dist)) - repmat(mean(dist_mtx(trl, 1:length(L_dist))), size(dist_mtx(trl, 1:length(L_dist))));
            %dist_mtx(trl, 1:length(L_dist)) = dist_mtx(trl, 1:length(L_dist))./std(dist_mtx(trl, 1:length(L_dist)));
        
            %plot distances
            %{
            hold on; 
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), R_dist, 'b--')
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), L_dist, 'g--')
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), dist_mtx(trl, 1:length(L_dist)), 'k')
            
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(R_dist, 50), 'b--')
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(L_dist, 50), 'g--')
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(dist_mtx(trl, 1:length(L_dist)), 50), 'k')
            title('left')
            
            %plot events
            plot([0 0], [-2 2], 'r', 'LineWidth', 3)  
            plot([stem_ent(trl) stem_ent(trl)], [-2 2], 'b', 'LineWidth', 3)
            plot([stem_ext(trl) stem_ext(trl)], [-2 2], 'b', 'LineWidth', 3)
            plot([rwd_time(trl) rwd_time(trl)], [-2 2], 'k', 'LineWidth', 3)
            %}
        
        %when right trial
        elseif trial_type(trl) == 2
            %higher values when right is smaller (inst rep is similar to right)
            dist_mtx(trl, 1:length(R_dist)) = L_dist - R_dist; 
            %dist_mtx(trl, 1:length(L_dist)) = dist_mtx(trl, 1:length(L_dist)) - repmat(mean(dist_mtx(trl, 1:length(L_dist))), size(dist_mtx(trl, 1:length(L_dist))));
            %dist_mtx(trl, 1:length(L_dist)) = dist_mtx(trl, 1:length(L_dist))./std(dist_mtx(trl, 1:length(L_dist)));
        
            %plot distances
            %{
            hold on; 
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), R_dist, 'b--')
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), L_dist, 'g--')
            plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), dist_mtx(trl, 1:length(L_dist)), 'k')
            
            
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(R_dist, 50), 'b--')
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(L_dist, 50), 'g--')
            %plot(time_set_align(trl, 1:sum(~isnan(time_set_align(trl, :)))), smooth(dist_mtx(trl, 1:length(L_dist)), 50), 'k')
            title('right')
            
            %plot events
            plot([0 0], [-2 2], 'r', 'LineWidth', 3)  
            plot([stem_ent(trl) stem_ent(trl)], [-2 2], 'b', 'LineWidth', 3)
            plot([stem_ext(trl) stem_ext(trl)], [-2 2], 'b', 'LineWidth', 3)
            plot([rwd_time(trl) rwd_time(trl)], [-2 2], 'k', 'LineWidth', 3)
            %}
            
        end
        %}
    end
    

    
    %light on bins
    align_pts = nan(size(time_set_align, 1), 1);
    for row = 1:size(time_set_align,1)     
        
        %skip bad trials and errors
        if ~ismember(row, incl_trials)
            continue
        end
        
        align_pts(row) = find(ismember(time_set_align(row, :), [-min(abs(time_set_align(row, :))) min(abs(time_set_align(row, :)))] ));
    end
  
    %align dist_mtx
    shifts = repmat(arb, size(align_pts)) - align_pts;
    for row = 1:size(dist_mtx,1)
        
        %skip bad trials and errors
        if ~ismember(row, incl_trials)
            continue
        end
        
        dist_mtx(row, :) = circshift(dist_mtx(row, :), shifts(row), 2);
    end
    
    %
    figure; imagesc(dist_mtx); 
    hold on
    plot([arb arb], [1 size(align_pts, 1)], 'r', 'LineWidth', .5)
    axis([find(nansum(dist_mtx)>0, 1, 'first') find(nansum(dist_mtx)>0, 1, 'last') .5 size(align_pts, 1)+.5])
    title('postshift')
   
    
    figure; plot(nanmean(dist_mtx)); 
    hold on
    plot([arb arb], [-2 2], 'r', 'LineWidth', 3)
    plot([find(nansum(dist_mtx)>0, 1, 'first') find(nansum(dist_mtx)>0, 1, 'last')], [0 0], 'k--', 'LineWidth', 1)
    axis([find(nansum(dist_mtx)>0, 1, 'first') find(nansum(dist_mtx)>0, 1, 'last') -2 2])
    %}
    

end