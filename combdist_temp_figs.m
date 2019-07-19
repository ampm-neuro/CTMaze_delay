
fig_name = comb_distmtx(nanmean(comb_distmtx_nonstd,2) > -10^14, :);

figure; imagesc(fig_name); 
    hold on
    plot([25000 25000], [.5 790.5], 'r', 'LineWidth', .5)
    axis([find(nansum(fig_name)>0, 1, 'first') find(nansum(fig_name)>0, 1, 'last') .5 size(fig_name,1)+.5])

    
    figure;

    xaxis = ((1:length(nanmean(fig_name)))-repmat(25000, size(1:length(nanmean(fig_name)))))./100;
    yaxis = smooth(nanmean(fig_name),20);
    
    apriori_correction = abs(mean(yaxis(xaxis> -5 & xaxis  < 0)));
    
plot(xaxis, yaxis+repmat(apriori_correction, size(yaxis)));
    hold on
    plot([0 0], [-2 2], 'r', 'LineWidth', 3)
    plot([-25 25], [0 0], 'k--', 'LineWidth', 1)
    axis([-5 5 -1 1])

    plot(xaxis, yaxis+smooth(nanstd(fig_name)./sqrt(sum(~isnan(fig_name))), 20)+repmat(apriori_correction, size(yaxis)), 'b-');
    plot(xaxis, yaxis-smooth(nanstd(fig_name)./sqrt(sum(~isnan(fig_name))), 20)+repmat(apriori_correction, size(yaxis)), 'b-');
    
    