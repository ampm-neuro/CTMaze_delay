function [summary interaction_only_count interaction_and_main_count main_only_count sect_dist_rank outcomes pvalue_lookup] = splitter_test(eptrials, clusters, lost_trials, out_stem, light_on_sect, accept_light_on, rat, session, idx)
%performs an emma wood style anova-based splitter test on each cell in clusters.
%
%output matrix: [cell_num, sig_sector, sig_trialtype, sig_interaction]

left_green=[52 153 70]./255;
right_blue=[46 49 146]./255;

sect_dist_rank = [];

%identify rejected trials
rejected_trials = reject_trials('stem', lost_trials, out_stem, light_on_sect, accept_light_on);

%trial choice types
%class_mtx = trl_class(eptrials, rejected_trials);
class_mtx = trl_class(eptrials, idx);

class_mtx = class_mtx(~isnan(class_mtx(:,3)),:);

%preallocate outcomes
outcomes = nan(length(clusters), 4);
pvalue_lookup = nan(length(clusters), 6);

%calculate firing rates eat each stem section
rate_mtx = stem_section_rates(eptrials, clusters, rejected_trials, idx);


%DOUBLE TRIAL TEST
%class_mtx = cat(1,class_mtx, class_mtx);
%rate_mtx = cat(1,rate_mtx, rate_mtx);


for clust = 1:length(clusters)
    cell = clusters(clust);
    
    left_rates = rate_mtx(class_mtx(:,3)==1 & class_mtx(:,4)==1, :, clust);
    right_rates = rate_mtx(class_mtx(:,3)==2 & class_mtx(:,4)==1, :, clust);
    
    %one repeated measures two-way anova, with sector as the rm
    p = anova_rm({left_rates right_rates}, 'off');
    outcomes(clust, :) = [cell p(1:3)<.05];
    
    %for figures
    pvalue_lookup(clust, :) = [rat, session, cell, p(1:3)];
    
    %1 if significant trialtype or trialtypeXsector interaction
    summary = sum(outcomes(:,3:4),2)>0;
    interaction_only_count = sum(outcomes(:,3:4),2)>0 & outcomes(:,4)==0;
    interaction_and_main_count = sum(outcomes(:,3:4),2)==2;
    main_only_count = sum(outcomes(:,3:4),2)>0 & outcomes(:,3)==0;
    
    %{
    figure
    hold on
    h_L = errorbar(mean(left_rates), std(left_rates)./sqrt(size(left_rates,1)), '-', 'linewidth', 2, 'Color', left_green); 
    h_R = errorbar(mean(right_rates), std(right_rates)./sqrt(size(right_rates,1)), '-', 'linewidth', 2, 'Color', right_blue); 
    legend([h_L h_R], 'Left Trials', 'Right Trials', 'location', 'northeastoutside');
    xlim([.5 4.5])
    set(gca, 'XTickLabel',{'Section 1','Section 2', 'Section 3', 'Section 4'},'XTick', 1:4, 'fontsize', 12)
    hold off
    yl = ylim;
    ylim([0 yl(2)])
    
    means = [mean(left_rates); mean(right_rates)];
    stds = [std(left_rates)./sqrt(size(left_rates,1)); std(right_rates)./sqrt(size(right_rates,1))];
    figure
    barweb(means(:), stds(:), .8)
    ylim([0 yl(2)])
    
   
    figure
    barweb(mean(left_rates), std(left_rates)./sqrt(size(left_rates,1)), .8)
    ylim([0 y1(2)])
    title('left')
    figure
    barweb(mean(right_rates), std(right_rates)./sqrt(size(right_rates,1)), .8)
    ylim([0 y1(2)])
    title('right')
    %}



end

%pos hoc ish
    interact_clusts = clusters(interaction_only_count | interaction_and_main_count);
    sect_dist_rank = nan(length(interact_clusts), 4);
    for clust = 1:length(interact_clusts)
        section_dists = nan(4,2);
        section_dists(:,1) = 1:4;
        for sect = 1:4
            section_dists(sect,2) = standard_distance(left_rates(:,sect), right_rates(:,sect));
        end
        
        section_dists = sortrows(section_dists, 2);
        sect_dist_rank(clust, :) = section_dists(:,1)';
    end
    
    
end