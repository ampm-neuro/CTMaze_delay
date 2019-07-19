function [mean_dif, mean_tstat] = stempos(eptrials, excluded, idx, cm_per_matlab_unit)
%calculates how long it took the rat to run down the stem on every trial

%trials to output
trials = setdiff(unique(eptrials(:,6))', excluded);
trials = trials(idx);

stem_pos = nan(length(trials), 2);


%iterate through trials
for i = 1:length(trials)
    
    trial = trials(i);
    
    stem_pos(i, 1:4) = [mean(eptrials(eptrials(:,6)==trial & eptrials(:,11)==1, 2))...
        mean(eptrials(eptrials(:,6)==trial & eptrials(:,11)==2, 2))...
        mean(eptrials(eptrials(:,6)==trial & eptrials(:,11)==3, 2))...
        mean(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4, 2))]; %pos in each sector
    stem_pos(i, 5) = mode(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11),1:4), 8)); %choice
    
end

    mean_dif_1 = abs(mean(stem_pos(stem_pos(:,5)==1, 1)) - mean(stem_pos(stem_pos(:,5)==2, 1)));
    mean_dif_2 = abs(mean(stem_pos(stem_pos(:,5)==1, 2)) - mean(stem_pos(stem_pos(:,5)==2, 2)));
    mean_dif_3 = abs(mean(stem_pos(stem_pos(:,5)==1, 3)) - mean(stem_pos(stem_pos(:,5)==2, 3)));
    mean_dif_4 = abs(mean(stem_pos(stem_pos(:,5)==1, 4)) - mean(stem_pos(stem_pos(:,5)==2, 4)));
    
    [~, ~, ~, stats1] = ttest2(stem_pos(stem_pos(:,5)==1, 1), stem_pos(stem_pos(:,5)==2, 1));
    [~, ~, ~, stats2] = ttest2(stem_pos(stem_pos(:,5)==1, 2), stem_pos(stem_pos(:,5)==2, 2));
    [~, ~, ~, stats3] = ttest2(stem_pos(stem_pos(:,5)==1, 3), stem_pos(stem_pos(:,5)==2, 3));
    [~, ~, ~, stats4] = ttest2(stem_pos(stem_pos(:,5)==1, 4), stem_pos(stem_pos(:,5)==2, 4));
    
    tstat1 = abs(stats1.tstat);
    tstat2 = abs(stats2.tstat);
    tstat3 = abs(stats3.tstat);
    tstat4 = abs(stats4.tstat);
    
    mean_dif = mean([mean_dif_1 mean_dif_2 mean_dif_3 mean_dif_4])*cm_per_matlab_unit;
    mean_tstat = mean([tstat1 tstat2 tstat3 tstat4]);
    

end