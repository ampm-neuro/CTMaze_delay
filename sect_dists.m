function [sds] = sect_dists(eptrials, clusters, lost_trials, light_on_sect, accept_lights)
%calculates for each cell in clusters the standard distance between L and 
%R trials at each corresponding maze section

%accuracy
acc = 1;

%accepted trials - nix lost trials and trials where the light doesnt turn on
rejected_trials = unique([lost_trials; light_on_sect(not(ismember(light_on_sect(:,2),accept_lights)),1)]);
%
%trial choice types
class_mtx = trl_class(eptrials, rejected_trials);

%number of sections
num_sects = 5;

%check number of trials, quit if too few
num_L_trials = length(class_mtx(class_mtx(:,3)==1 & class_mtx(:,4)==acc, 1));
num_R_trials = length(class_mtx(class_mtx(:,3)==2 & class_mtx(:,4)==acc, 1));
if num_L_trials<3 || num_R_trials<3
    sds = nan(length(clusters), num_sects);
    display('Too few of one or both trial types. Terminating sect_dists.')
    return
end
    
%rate matrix
[~, msrs_out] = maze_section_rates(eptrials, clusters, lost_trials, light_on_sect, [1 2.1]);
if size(msrs_out,2) == 8
    msrs = msrs_out(:, 1:5, :);
end

%msrs

%preallocate outcomes
sds = nan(length(clusters), num_sects);


if size(msrs, 1) ~= size(class_mtx,1)
    num_trials_rate_mtx = size(msrs, 1)
    num_trials_class_mtx = size(class_mtx,1)
    error('rate matrix (msrs) should have the same number of rows (trials) as the class matrix (class_mtx)')
end


%find mean and pooled sd from msrs
    for row = 1:length(clusters)        
        for sect = 1:num_sects
            
            %calculate pooled standard deviation
            %
                %section rates
                
                if size(msrs,2)==7
                    if ismember(sect,1:3)
                        r_sL = msrs(class_mtx(:,3)==1 & class_mtx(:,4)==acc, sect, row);
                        r_sR = msrs(class_mtx(:,3)==2 & class_mtx(:,4)==acc, sect, row);
                    elseif sect==4
                        r_sL = msrs(class_mtx(:,3)==1 & class_mtx(:,4)==acc, 4, row);
                        r_sR = msrs(class_mtx(:,3)==2 & class_mtx(:,4)==acc, 5, row);
                    elseif sect==5
                        r_sL = msrs(class_mtx(:,3)==1 & class_mtx(:,4)==acc, 6, row);
                        r_sR = msrs(class_mtx(:,3)==2 & class_mtx(:,4)==acc, 7, row);
                    end
                elseif size(msrs,2)==5
                    r_sL = msrs(class_mtx(:,3)==1 & class_mtx(:,4)==acc, sect, row);
                    r_sR = msrs(class_mtx(:,3)==2 & class_mtx(:,4)==acc, sect, row);
                else
                    size(msrs)
                    error('check msrs size')
                end
                
                %section means
                m_sL = mean(r_sL);
                m_sR = mean(r_sR);

                %section variances
                var_sL = var(r_sL);
                var_sR = var(r_sR);

                %pooled var and std
                pooled_var = sum([var_sL*(length(r_sL-1)) var_sR*(length(r_sR-1))])/(length([r_sL; r_sR])-2);
                %pooled_var = sum([var_sL*(length(r_sL)) var_sR*(length(r_sR))])/(length([r_sL; r_sR]));
                pooled_std = sqrt(pooled_var);
                
            %load standard difference
            sds(row, sect) = abs(m_sL-m_sR)/pooled_std;
            
        end
    end
end

