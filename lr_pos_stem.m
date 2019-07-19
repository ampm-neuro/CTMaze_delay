function lr_pos_stem(eptrials, section_boundaries)
%plot stem trajectories in green for left choice and blue for right choice

grn=[52 153 70]./255;
blu=[46 49 146]./255;

figure
hold on
for trial= 1:max(eptrials(:,6))
    if mode(eptrials(eptrials(:,6)==trial,8))==1
    
        plot(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11), 1:4),2), eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11), 1:4),3), 'Color', grn, 'LineWidth', 0.5, 'LineStyle', '-');
        
    elseif mode(eptrials(eptrials(:,6)==trial,8))==2
        
        plot(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11), 1:4),2), eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11), 1:4),3), 'Color', blu, 'LineWidth', 0.5, 'LineStyle', '-');

    end
end


sections(eptrials, section_boundaries)


end