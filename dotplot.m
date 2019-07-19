function dotplot(eptrials, cluster)
%makes a dot plot for left and right trials, side by side. spacing between
%them determined by the variable 'spacing'. try 100.

spacing = 100;


figure
set(gca,'xdir','reverse')
hold on


%trajectories
for trial= 1:max(eptrials(:,6))
    
    
    local_eptrials = eptrials(eptrials(:,6)==trial,:);
    last_iti = max(local_eptrials(local_eptrials(:,10)==8, 1));
    local_eptrials(ismember(local_eptrials(:,10),[8 9]), :) = nan;
    if ~isempty(last_iti)
        local_eptrials(local_eptrials(:,1)>last_iti, :) = nan;
    end
    
    
    if mode(local_eptrials(:,8))==1
        plot(local_eptrials(:,2), local_eptrials(:,3), 'Color', [.8 .8 .8], 'LineWidth', 0.5, 'LineStyle', '-');
        
    elseif mode(local_eptrials(:,8))==2
        plot(local_eptrials(:,2) + repmat(-spacing, size(local_eptrials(:,2))), local_eptrials(:,3), 'Color', [.8 .8 .8], 'LineWidth', 0.5, 'LineStyle', '-');

    end
end

%spikes
for trial= 1:max(eptrials(:,6))
    
    
    local_eptrials = eptrials(eptrials(:,6)==trial,:);
    last_iti = max(local_eptrials(local_eptrials(:,10)==8, 1));
    local_eptrials(ismember(local_eptrials(:,10),[8 9]), :) = nan;
    if ~isempty(last_iti)
        local_eptrials(local_eptrials(:,1)>last_iti, :) = nan;
    end
    
    
    if mode(local_eptrials(:,8))==1
        plot(local_eptrials(local_eptrials(:,4)==cluster,2), local_eptrials(local_eptrials(:,4)==cluster,3), 'r.', 'Markersize', 6);
        
    elseif mode(local_eptrials(:,8))==2
        plot(local_eptrials(local_eptrials(:,4)==cluster,2) + repmat(-spacing, size(local_eptrials(local_eptrials(:,4)==cluster,2))), local_eptrials(local_eptrials(:,4)==cluster,3), 'r.', 'Markersize', 6);

    end
end

xlim([700 1200])
axis off

end