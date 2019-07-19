function lr_pos(eptrials, section_boundaries)
%plot non-carry trajectories in green for left choice and blue for right choice

grn=[52 153 70]./255;
blu=[46 49 146]./255;

figure
hold on
for trial= 1:max(eptrials(:,6))
    
    
    local_eptrials = eptrials(eptrials(:,6)==trial,:);
    last_iti = max(local_eptrials(local_eptrials(:,10)==8, 1));
    local_eptrials(local_eptrials(:,10)==9, 2:3) = nan;
    if ~isempty(last_iti)
        local_eptrials(local_eptrials(:,1)>last_iti, 2:3) = nan;
    end
    
    
    if mode(local_eptrials(:,8))==1
    
        plot(local_eptrials(:,2), local_eptrials(:,3), 'Color', grn, 'LineWidth', 0.5, 'LineStyle', '-');
        
    elseif mode(local_eptrials(:,8))==2
        
        plot(local_eptrials(:,2), local_eptrials(:,3), 'Color', blu, 'LineWidth', 0.5, 'LineStyle', '-');
    end
end


sections(eptrials, section_boundaries)


end