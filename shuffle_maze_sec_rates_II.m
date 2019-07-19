function [shuf_dists shuf_stds shuf_sems shuf_mean_sems] = shuffle_maze_sec_rates_II(number_of_shuffles, all_rates, all_LRs)
%takes output from ALL_maze_sec_rates, which is in the form of 3 cell
%arrays, each containing 2 cell arrays, one of which contains 1 3d matrix
%per session and the other contains a list of trial choices for each
%session
%
%shuffles the trial choices in each choice session cell, or not
%
%uses trial choices in each choice session cell to index each 3d session
%matrix, and calculate difference between L and R means at each section
%(possibly divided by the pooled sd, or something)
%
%outputs the average distance at each section over all cells

% stage number
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day

%zero number of shuffles means use actual data
if number_of_shuffles == 0
    num_shuffles = 1;
    shuf_stds = nan(4,5);
    shuf_sems = nan(4,5);
    shuf_mean_sems = nan(4,1);
else
    num_shuffles = number_of_shuffles;
    shuf_stds = [];
    shuf_sems = [];
    shuf_mean_sems = [];
end



%TEST CONSTRICT TRIAL NUMBERS TO 20
%{
for istg = 1:size(all_rates,2)
    for isesh = 1:size(all_LRs{istg},2)
            all_rates{istg}{isesh} = all_rates{istg}{isesh}(end-19:end, :,:);
            all_LRs{istg}{isesh} = all_LRs{istg}{isesh}(end-19:end);
    end
end
%}            

%preallocate output
dists = nan(size(all_rates,2),5);
concat_dists = cell(1,size(all_rates,2));
mean_dists = cell(1,size(all_rates,2));
shuf_dists{1} = nan(num_shuffles, 5);
shuf_dists{2} = nan(num_shuffles, 5);
shuf_dists{3} = nan(num_shuffles, 5);

%iterate through shuffles at least once (even if none requested)
for shuf = 1:num_shuffles

    %shuffle LR's
    shuffled_LRs = cell(1,size(all_rates,2));
    for stge_shuf = 1:size(all_rates,2)
        
        shuffled_LRs{stge_shuf} = cell(1,size(all_LRs{stge_shuf},2));
        for sesh_shuf = 1:size(all_LRs{stge_shuf},2)
            
            %fill a mirror of all_LRs with the corresponding session
            %indexed by a random permutation
            shuffled_LRs{stge_shuf}{sesh_shuf} = all_LRs{stge_shuf}{sesh_shuf}(randperm(length(all_LRs{stge_shuf}{sesh_shuf})));

        end
    end

    %calculate rates
    for stage = 1:size(all_rates,2)

        %current_stage = stage        
        for session = 1:size(all_rates{stage},2)
            
            %current_session = session
            %index out rates from left and right (correct only) trials
            if number_of_shuffles == 0 %keep original LRs if no shuffles
                left_rates = all_rates{stage}{session}(all_LRs{stage}{session}==1, :, :);
                right_rates = all_rates{stage}{session}(all_LRs{stage}{session}==2, :, :);
            else %otherwise use shuffled LRs
                left_rates = all_rates{stage}{session}(shuffled_LRs{stage}{session}==1, :, :);
                right_rates = all_rates{stage}{session}(shuffled_LRs{stage}{session}==2, :, :);
            end

            mean_left = nanmean(left_rates);
            mean_right = nanmean(right_rates);

            std_left = nanstd(left_rates);
            std_right = nanstd(right_rates);

            %can change this method
            std_crude_pool = (std_left+std_right)./2;
            %std_crude_pool = std_crude_pool./std_crude_pool;%mean only
            %std_crude_pool = std_left*(size(left_rates,3)/(size(left_rates,3)+size(right_rates,3))) + std_right*(size(right_rates,3)/(size(left_rates,3)+size(right_rates,3)));%weighted average std

            mean_differences = abs(mean_left - mean_right);

            std_dists = mean_differences./std_crude_pool;

            %load stage cell
            concat_dists{stage} = [concat_dists{stage}; reshape(squeeze(std_dists)', size(std_dists,3), size(std_dists,2))];
            %concat_dists{stage} = [concat_dists{stage};%nanmean(reshape(squeeze(std_dists)', size(std_dists,3),%size(std_dists,2)),1)];%avg within session first

            %load means
            mean_dists{stage} = [mean_dists{stage}; nanmean(reshape(squeeze(std_dists)', size(std_dists,3), size(std_dists,2)),2)];
        end
        
        %load stage mean into dists output
        dists(stage,:) = nanmean(concat_dists{stage});
        
        
        if number_of_shuffles == 0            
            shuf_stds(stage,:) = nanstd(concat_dists{stage});
            shuf_sems(stage,:) = nanstd(concat_dists{stage})./sqrt(size(concat_dists{stage},1));
            shuf_mean_sems(stage,:) = nanstd(mean_dists{stage})./sqrt(size(mean_dists{stage},1));
        end
        
    end
    
    %load current shuffle outcomes
    for stg = 1:size(all_rates,2)
        shuf_dists{stg}(shuf, :) = dists(stg,:);
    end
    
end