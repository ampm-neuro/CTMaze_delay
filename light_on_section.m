function light_on_sect = light_on_section(eptrials, section_boundaries, lost_trials)
%determines the rats location when the light came on on each trial
%
%light_on_sect indicates the section the rat was in when the light came on
%
% 0 = light did not stay on
% 1 = start area
% 2.1 = first stem section
% 2.2 = second stem section
% 2.3 = third stem section
% 2.4 = fourth stem section
% 3 = choice area
% 4 = left arm
% 5 = right arm
% 6 = left reward section
% 7 = right reward section
% 8 = iti section
%


%stem boundaries
    %center of maze
    comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
    comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;
    %stem section boundaries [xlow xhigh ylow yhigh]
    stem(1,:) = [comx-30 comx+30 comy-90 comy-52.5];
    stem(2,:) = [comx-30 comx+30 comy-52.5 comy-15];
    stem(3,:) = [comx-30 comx+30 comy-15 comy+22.5];
    stem(4,:) = [comx-30 comx+30 comy+22.5 comy+60];

%preallocate
light_on_sect = zeros(max(eptrials(:,6)),1);

%iterate through trials
for trial = setdiff(1:max(eptrials(:,6)), lost_trials)
    
    %if the light doesnt stay on, just leave it as 0
    if length(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13), 1:2), 1)) > 10 %probably just needs to be > 1
    
        light_on = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13), 1:2), 1));

        rat_pos = eptrials(eptrials(:,1)==light_on, 2:3);

        %match rat_pos with a row from section_boundaries
        for row = 1:size(section_boundaries(:,1))
            if rat_pos(1) > section_boundaries(row,1) && rat_pos(1) < section_boundaries(row,2) && rat_pos(2) > section_boundaries(row,3) && rat_pos(2) < section_boundaries(row,4)
                light_on_sect(trial) = row;
                break
            end
        end

        %if it came on in stem, indicate which stem section using decimal
        %notation
        if light_on_sect(trial) == 2
            for row = 1:size(section_boundaries(:,1))
                if rat_pos(1) > stem(row,1) && rat_pos(1) < stem(row,2) && rat_pos(2) > stem(row,3) && rat_pos(2) < stem(row,4)
                    light_on_sect(trial) = light_on_sect(trial) + row/10;
                    break
                end
            end
        end  
    end
end

light_on_sect = [(1:max(eptrials(:,6)))' light_on_sect];

end