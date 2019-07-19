function rewards(eptrials)
%plot reward flag locations over grey trajectory

figure
plot(eptrials(:,2), eptrials(:,3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-');
hold on

rwd_times = eptrials(ismember(eptrials(:,5), [1 2]), 1);

for rwd = 1:length(rwd_times)
    
    if eptrials(eptrials(:,1)==rwd_times(rwd), 5) == 1
        
        plot(eptrials(eptrials(:,1)==rwd_times(rwd), 2), eptrials(eptrials(:,1)==rwd_times(rwd), 3), 'b.', 'markersize', 10); 
        
    elseif eptrials(eptrials(:,1)==rwd_times(rwd), 5) == 2
        
        plot(eptrials(eptrials(:,1)==rwd_times(rwd), 2), eptrials(eptrials(:,1)==rwd_times(rwd), 3), 'r.', 'markersize', 10);
        
    end
end

end