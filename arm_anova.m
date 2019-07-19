function [summary pvalue_lookup] = arm_anova(eptrials, clusters, rat, session)
%look for differential finding on the non-overlapping sections (choice,
%arm, rwd)

%preal
summary = nan(length(clusters), 2);
pvalue_lookup = nan(length(clusters), 5);
count = 0;

%iterate through clusters
for clust = clusters'
    
    count = count + 1;
    
    %choice, arm, rwd section rates
    rate_mtx = maze_sec_rates_II(eptrials, clust);
    left_rates = rate_mtx(rate_mtx(:,6)==1 & rate_mtx(:,7)==1,3:6);
    right_rates = rate_mtx(rate_mtx(:,6)==2 & rate_mtx(:,7)==1,3:6);

    
    left_choice = left_rates(:,1);
    left_arm = left_rates(:,2);
    right_choice = right_rates(:,1);
    right_arm = right_rates(:,2);
    
    [h_choice p_choice] = ttest2(left_choice, right_choice);
    [h_arm p_arm] = ttest2(left_arm, right_arm);
    
    summary(count, :) = [h_choice h_arm];
    
    
    %anova input vectors
    %{
    rates = [left_rates(:,1:3);right_rates(:,1:3)];
    type = [left_rates(:,4) left_rates(:,4) left_rates(:,4);right_rates(:,4) right_rates(:,4) right_rates(:,4)];
    section = repmat(1:3, size(rates,1), 1);
    
    
    %anova
    p = anovan(rates(:), {type(:) section(:)}, 'model', 2, 'sstype', 2, 'varnames', char('LR', 'sect'), 'display', 'off')';
    
    outcomes(count, :) = [clust p(1:3)<.05];
    %}

    
    cur_var_L = left_arm;
    cur_var_R = right_arm;
    
    
    left_green=[52 153 70]./255;
    right_blue=[46 49 146]./255;
    figure
    hold on
    h_L = errorbar(mean(cur_var_L), std(cur_var_L)/sqrt(length(cur_var_L)), '-', 'linewidth', 2, 'Color', left_green);
    h_R = errorbar(mean(cur_var_R), std(cur_var_R)/sqrt(length(cur_var_R)), '-', 'linewidth', 2, 'Color', right_blue); 
    legend([h_L h_R], 'Left Trials', 'Right Trials', 'location', 'northeastoutside');
    xlim([.95 1.05])
    yL = ylim;
    ylim([0 yL(2)])

    hold off
    figure
    barweb([mean(cur_var_L) mean(cur_var_R)], [std(cur_var_L)/sqrt(length(cur_var_L)) std(cur_var_R)/sqrt(length(cur_var_R))], .8)
    ylim([0 yL(2)])
    
    pvalue_lookup(count, :) = [rat, session, clust, p_choice, p_arm];

end

%summary = sum(outcomes(:,3:4),2)>0;