function [summary outcomes pvals] = light_rwd_cell(eptrials, clusters, windowback, windowfwd, plot_yn)
%calculates firing rates aroudn the light-on and reward on each trial.
%performs an 2x2 trial type by before/after anova on each the lighton and
%the reward, and reports signif or no in summary and whether it was an
%interaction or w/e in outcomes.

%preallocate
trl_rates = nan(max(eptrials(:,6)), 4, length(clusters));
outcomes = nan(length(clusters), 6);
pvals = nan(size(outcomes));
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
            first_rwd = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,5), 1:4) ,1));
            FR_time(trial) = first_rwd;
        end

        %firing rates
        time_before_light = eptrials(:,1) > FL_time(trial)-windowback & eptrials(:,1) < FL_time(trial);
        spikes_before_light = sum(eptrials(:,4)==cluster & time_before_light);
        FR_before_light = spikes_before_light/windowback;

        time_after_light = eptrials(:,1) > FL_time(trial) & eptrials(:,1) < FL_time(trial)+windowfwd;
        spikes_after_light = sum(eptrials(:,4)==cluster & time_after_light);
        FR_after_light = spikes_after_light/windowfwd;

        time_before_rwd = eptrials(:,1) > FR_time(trial)-windowback & eptrials(:,1) < FR_time(trial);
        spikes_before_rwd = sum(eptrials(:,4)==cluster & time_before_rwd);
        FR_before_rwd = spikes_before_rwd/windowback;

        time_after_rwd = eptrials(:,1) > FR_time(trial) & eptrials(:,1) < FR_time(trial)+windowfwd;
        spikes_after_rwd = sum(eptrials(:,4)==cluster & time_after_rwd);
        FR_after_rwd = spikes_after_rwd/windowfwd;


        %load trial rates
        trl_rates(trial,:, row) = [FR_before_light FR_after_light FR_before_rwd FR_after_rwd];

        
    end

    %find trial accuracy
    class_idx = trl_class(eptrials);
    
    if row==1
    trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1,1:2,1);
    end


    %perform anovas; within subjects time, between subjects trialtype
    %output pvalues in order: before/after, trialtype, interaction
        %light
        p_light = anova_rm({[trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 2, row)] [trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 1, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 2, row)]}, 'off');
        %rwd
        p_rwd = anova_rm({[trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==1, 4, row)] [trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 3, row) trl_rates(class_idx(:,4)==1 & class_idx(:,3)==2, 4, row)]}, 'off');

        %light (before/after, trialtype, interaction) rwd (before/after,
        %trialtype, interaction)
        outcomes(row, :) = [p_light(1:3)<.05 p_rwd(1:3)<.05];
        pvals(row, :) = [p_light(1:3) p_rwd(1:3)];

        %light cell if interaction or main effect time; discriminating light cell
        %if interaction or.. main effect time AND main effect trialtype
        %
        %rwd cell if interaction or main effect time; discriminating rwd cell
        %if interaction or main effect time AND main effect trialtype
        %
        %[light_cell light_discrim_cell rwd_cell reward_discrim_cell]
        summary(row, :) = [[sum(outcomes(row,[1 3]))>0 sum([outcomes(row,3)==1 sum(outcomes(row,1:2))==2]) > 0] [sum(outcomes(row,[4 6]))>0 sum([outcomes(row,6)==1 sum(outcomes(row,4:5))==2]) > 0]];

    

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