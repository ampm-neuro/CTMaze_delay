function newvid = iso_pos_scrub(vid)

%add time and position (columns 1:3)
    epvid = vid;
    a = 1:length(epvid);
%test figures
    figure; plot(epvid(:,2), epvid(:,3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-')
    title('original')
    figure; plot3(epvid(:,2), epvid(:,3), a, 'Color', [0.9 0.9 0.9] , 'LineWidth', 0.5, 'LineStyle', '-')
    hold on
    plot3(epvid(:,2), epvid(:,3), a, '.', 'Color', [0.8 0.8 0.8], 'markersize', 6)

   %remove improbable changes in position
    too_big_origin = 250;
    too_big = too_big_origin;
    cum_mod = 5;
    
    p1 = epvid(1, 2:3);
    p2 = epvid(2, 2:3);
  
    time_standard = epvid(:,1)./1000000;
    
    t1 = time_standard(1,1);
    t2 = time_standard(2,1);
    
    deletions = zeros(length(epvid(:,2)),1);
    dists = zeros(length(epvid(:,2)),1);

    count = 0;
    for i = 1:length(epvid(:,2))-2
        
        current_distance = pdist([p1; p2])/(t2-t1);
        dists(i+1) = current_distance;
        
        if current_distance > too_big
            deletions(i+1) = 1;
            
            p2 = epvid(i+2, 2:3);
            t2 = time_standard(i+2, 1);
            
            count = count + cum_mod;
            too_big = too_big + count;
        
        else
            too_big = too_big_origin;
            count = 0;
            
            p1 = epvid(i+1, 2:3);
            p2 = epvid(i+2, 2:3);
            
            t1 = time_standard(i+1, 1);
            t2 = time_standard(i+2, 1);
            
        end
    end
    
    deletions = logical(deletions);
    
    
    plot3(epvid(deletions,2), epvid(deletions,3), a(deletions), '.', 'Color', [1 0 0], 'markersize', 6)
    title('dubious points')
    
    epvid(deletions, 2:3) = NaN;
        
    figure; plot(epvid(~isnan(epvid(:,2)),2), epvid(~isnan(epvid(:,2)),3), 'Color', [0.8 0.8 0.8] , 'LineWidth', 0.5, 'LineStyle', '-')
    title('after deletion')
    
    newvid = epvid;
    

end