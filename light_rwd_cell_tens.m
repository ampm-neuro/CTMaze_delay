function [outcomes_2groups] = light_rwd_cell_tens(eptrials, clusters, windowback, windowfwd, plot_yn)
%calculates firing rates aroudn the light-on and reward on each trial.
%performs an 2x2 trial type by before/after anova on each the lighton and
%the reward, and reports signif or no in summary and whether it was an
%interaction or w/e in outcomes.

%preallocate
trl_rates = nan(max(eptrials(:,6)), 40, length(clusters));
outcomes_2groups = nan(length(clusters), 6);
outcomes_3groups = nan(length(clusters), 14);
summary = nan(length(clusters), 4);

%first light times
FL_time = nan(max(eptrials(:,6)), 1);

%first rwd times
FR_time = nan(max(eptrials(:,6)), 1);

for row = 1:length(clusters)
    cluster = clusters(row);

    %calculate rates
    for trial = 1:max(eptrials(:,6))

        if row == 1
            %first instant of light-on
            first_light = min(eptrials(eptrials(:,6)==trial & ~isnan(eptrials(:,13)),1));
            
            %if no light, use stem entrance (should only be 1 trial)
            if isempty(first_light)
                FL_time(trial) = min(eptrials(eptrials(:,6)==trial & eptrials(:,10)==2,1));
            else
                FL_time(trial) = first_light;
            end

            %first reward flag
            first_rwd = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,5), 1:2) ,1));
            FR_time(trial) = first_rwd;
        end

        FR_before_light = nan(10,1);
        FR_after_light = nan(10,1);
        FR_before_rwd = nan(10,1);
        FR_after_rwd = nan(10,1);
        
        
        for bins = 1:10
        
            %time bounds before
            current_bin_low = (bins-1)*(windowback/10);
            current_bin_hi = (bins)*(windowback/10);
            
            %firing rates before
            before_light = FL_time(trial)-windowback;
            wndw_low = before_light + current_bin_low;
            wndw_hi = before_light + current_bin_hi;         
            current_light_bin_time = eptrials(:,1) >= wndw_low & eptrials(:,1) < wndw_hi;
            spikes_in_current_light_bin = sum(eptrials(:,4)==cluster & current_light_bin_time);
            FR_before_light(bins,1) = spikes_in_current_light_bin/(wndw_hi - wndw_low);

            before_rwd = FR_time(trial)-windowback;
            wndw_low = before_rwd + current_bin_low;
            wndw_hi = before_rwd + current_bin_hi;
            current_light_bin_time = eptrials(:,1) >= wndw_low & eptrials(:,1) < wndw_hi;
            spikes_in_current_rwd_bin = sum(eptrials(:,4)==cluster & current_light_bin_time);
            FR_before_rwd(bins,1) = spikes_in_current_rwd_bin/(wndw_hi - wndw_low);
            
            
            %time bounds after
            current_bin_low = (bins-1)*(windowfwd/10);
            current_bin_hi = (bins)*(windowfwd/10);
            
            %firing rates after
            light_time = FL_time(trial);
            wndw_low = light_time + current_bin_low;
            wndw_hi = light_time + current_bin_hi;
            current_light_bin_time = eptrials(:,1) > wndw_low & eptrials(:,1) <= wndw_hi;
            spikes_in_current_light_bin = sum(eptrials(:,4)==cluster & current_light_bin_time);
            FR_after_light(bins,1) = spikes_in_current_light_bin/(wndw_hi - wndw_low);

            rwd_time = FR_time(trial);
            wndw_low = rwd_time + current_bin_low;
            wndw_hi = rwd_time + current_bin_hi;
            current_rwd_bin_time = eptrials(:,1) > wndw_low & eptrials(:,1) <= wndw_hi;
            spikes_in_current_rwd_bin = sum(eptrials(:,4)==cluster & current_rwd_bin_time);
            FR_after_rwd(bins,1) = spikes_in_current_rwd_bin/(wndw_hi - wndw_low);
            
            
        
        end
        
        %load trial rates 40 across e.g., trl_rates(:,1:10) = rates before
        %the light
        trl_rates(trial,:, row) = [FR_before_light' FR_after_light' FR_before_rwd' FR_after_rwd'];

    end
    
    %find trial accuracy
    class_idx = trl_class(eptrials);

    %perform anovas; within subjects time, between subjects trialtype
    %output pvalues in order: before/after, trialtype, interaction
        
    
    
    
        %light 2 groups (bins x trialtype)
        left_before = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 1:10, row);
        left_after = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 11:20, row);
        right_before = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 1:10, row);
        right_after = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 11:20, row);
        
        
        if row==1
        [left_before left_after]
        [mean(left_before,2) mean(left_after,2)]
        end

        p_light_2groups = anova_rm({[left_before left_after] [right_before right_after]}, 'off');
        
        
        %rwd 2 groups (bins x trialtype)
        left_before = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 21:30, row);
        left_after = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 31:40, row);
        right_before = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 21:30, row);
        right_after = trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 31:40, row);
        
        p_rwd_2groups = anova_rm({[left_before left_after] [right_before right_after]}, 'off');

        %light (before/after, trialtype, interaction) rwd (before/after,
        %trialtype, interaction)
        outcomes_2groups(row, :) = [p_light_2groups(1:3) p_rwd_2groups(1:3)];

        
        %{
        %light 3 groups (bins x before/after x trialtype)
        accepted_trials = class_idx(:,4)==1;
        light_data = trl_rates(accepted_trials, 1:20, row);
        subject = repmat((1:max(eptrials(:,6)))', 1, size(light_data,2));
        subject = subject(accepted_trials,:);
        subject = repmat((1:size(subject,1))', 1, size(light_data,2))
        bin = repmat([1:10 1:10], max(eptrials(:,6)), 1);
        bin = bin(accepted_trials,:)
        bef_aft = repmat([ones(1,10) repmat(2, 1, 10)], max(eptrials(:,6)), 1);
        bef_aft = bef_aft(accepted_trials,:)
        trl_type = repmat(class_idx(:,3), 1, size(light_data,2));
        trl_type = trl_type(accepted_trials,:)
        sum(trl_type==1)
        sum(trl_type==2)

        p_light_3groups = anova_rm_3([light_data(:) bef_aft(:) bin(:) trl_type(:) subject(:)], 0.05);
        
        
        %rwd
        %rwd_data = trl_rates(accepted_trials, 21:40, row);
        
        %p_rwd_3groups = anova_rm_3([rwd_data(:) bin(:) bef_aft(:) trl_type(:) subject(:)], 0.05);
        
        
        %light (before/after, trialtype, interaction) rwd (before/after,
        %trialtype, interaction)
        %outcomes_3groups(row, :) = [p_light_3groups(1:7)<.05 p_rwd_3groups(1:7)<.05];        
        %}
        
        
        
        
        %light cell if interaction or main effect time; discriminating light cell
        %if interaction or.. main effect time AND main effect trialtype
        %
        %rwd cell if interaction or main effect time; discriminating rwd cell
        %if interaction or main effect time AND main effect trialtype
        %
        %[light_cell light_discrim_cell rwd_cell reward_discrim_cell]
        %summary(row, :) = [[sum(outcomes(row,[1 3]))>0 sum([outcomes(row,3)==1 sum(outcomes(row,1:2))==2]) > 0] [sum(outcomes(row,[4 6]))>0 sum([outcomes(row,6)==1 sum(outcomes(row,4:5))==2]) > 0]];

    

%plot
if plot_yn ==1
    
    grn=[52 153 70]./255;
    blu=[46 49 146]./255;
    
    if row == 1
        %light locations
        figure
        hold on
        plot(eptrials(:, 2), eptrials(:, 3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-') 
        for trial = 1:max(eptrials(:,6))

            time_idx = eptrials(:,1)>FL_time(trial)-windowback & eptrials(:,1) < FL_time(trial)+windowfwd;

            if class_idx(trial,3)==1 && class_idx(trial,4)==1
                plot(eptrials(time_idx,2), eptrials(time_idx,3), 'Color', grn, 'LineWidth', 0.5, 'LineStyle', '-');
            elseif class_idx(trial,3)==2 && class_idx(trial,4)==1
                plot(eptrials(time_idx,2), eptrials(time_idx,3), 'Color', blu, 'LineWidth', 0.5, 'LineStyle', '-');
            end
        end
        hold off
        title('Light pos')

        %rwd locations
        figure
        hold on
        plot(eptrials(:, 2), eptrials(:, 3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-') 
        for trial = 1:max(eptrials(:,6))

            time_idx = eptrials(:,1)>FR_time(trial)-windowback & eptrials(:,1) < FR_time(trial)+windowfwd;

            if class_idx(trial,3)==1 && class_idx(trial,4)==1
                plot(eptrials(time_idx,2), eptrials(time_idx,3), 'Color', grn, 'LineWidth', 0.5, 'LineStyle', '-');
            elseif class_idx(trial,3)==2 && class_idx(trial,4)==1
                plot(eptrials(time_idx,2), eptrials(time_idx,3), 'Color', blu, 'LineWidth', 0.5, 'LineStyle', '-');
            end
        end
        hold off
        title('Rwd pos')
    end
    
    %light
    figure
    barvalues_L = mean([trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 2, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 2, row)]);
    errors_L = std([trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 2, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 2, row)])./sqrt(length(trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 2, row)));
    barweb(barvalues_L, errors_L, .8);
    title(strcat('Light ', num2str(cluster)))
    
    %rwd
    figure
    barvalues_R = mean([trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 4, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 4, row)]);
    errors_R = std([trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 4, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 4, row)])./sqrt(length(trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 4, row)));
    barweb(barvalues_R, errors_R, .8);
    title(strcat('Reward ', num2str(cluster)))
    
    
    
end
end


end