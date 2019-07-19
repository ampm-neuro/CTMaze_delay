function [eptrials section_boundaries lost_trials] = lv_pretrials(vid, ttun, spiketimes, REWARD_L, REWARD_R, ITI_ON, ITI_OFF, LIGHT_ON, FAKE_ON, LIGHT_OFF, TRIAL_START)
%   builds the super matrix 'eptrials' from lindsey's data matrices
%
% INPUT VARIABLES (18 total):
%   vid = timestamps and x and y positions
%   ttun = tetrode and unit number in rows
%   spiketimes = the spike timestamps for the tetrode and unit from ttun
%   REWARD_ L/R = first timestamp of licks on left/right error trials
%   ITI_ ON/OFF = first/last timestamp of intertrial interval
%   LIGHT_ ON/OFF = first/last timestamp of left light being on
%   TRIAL_START = first timestamp of every trial
%
% OUTPUT VARIABLE eptrials:
%   eptrials(:, 1) time
%   eptrials(:, 2) X position
%   eptrials(:, 3) Y position
%   eptrials(:, 4) is a spike (tetrode.cluster), video sample (1), or NaN
%   eptrials(:, 5) manual flag (rwd_L/R = 1/2)
%   eptrials(:, 6) trial number
%   eptrials(:, 7) NAN
%   eptrials(:, 8) choice (1 = L, 2 = R)
%   eptrials(:, 9) ones
%   eptrials(:,10) maze section
%       1 = start area
%       2 = stem
%       3 = choice area
%       4 = choice arm L
%       5 = choice arm R
%       6 = reward area L
%       7 = reward area R
%       8 = ITI platform
%       9 = carry
%   eptrials(:,11) stem section 1:4, ordered from start to choice
%   eptrials(:,12) velocity (cm/s)
%   eptrials(:,13) light status (1 = L light on, 2 = R light on, or nan)
%


%find neuralynx start time
initial_clock = min(vid(:,1));

%identify the trial numbers
trial_numbers = 1:length(TRIAL_START);

%initialize lost trials
lost_trials = [];

%create eptrials from time and position data (columns 1:3)
eptrials = time_and_pos(vid, 300, 5, initial_clock);
    function eptrials = time_and_pos(vid, too_big_origin, cum_mod, initial_clock)
    % time_and_pos attempts to correct position data by removing spurious
    % points creates by light interference. It also resamples the data at
    % 100hz.
    %
    % INPUT VARIABLES
    %   vid = the nueralynx video data with columns time, xpos, and ypos
    %   too_big_origin = the max acceptable position change between adjacent
    %      %points. Larger position changes are deleted as noise.
    %   cum_mod = a correction factor that prevents deletions from creating
    %      %unacceptable position changes between points
    %   initial_clock = neuralynx start time
    %
    % OUTPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   initial_clock = the uncorrected session start time
    %
    
    
        %create matrix eptrials from vid input
        eptrials = vid;

        %set start time to 0sec
        eptrials(:,1) = correct_time(eptrials(:,1), initial_clock);

        %prep to remove improbable changes in position
        too_big = too_big_origin;
        guilty_by_association = 4;

        %adjacent pos points for evaluating velocity
        p1 = eptrials(1, 2:3);
        p2 = eptrials(2, 2:3);

        %adjacent time points for evaluating velocity
        t1 = eptrials(1,1);
        t2 = eptrials(2,1);

        %preallocate
        deletions = zeros(length(eptrials(:,2)),1);
        dists = zeros(length(eptrials(:,2)),1);

        %iterate through adjacent points to evaluate velocity
        count = 0;
        for i = 1:length(eptrials(:,2))-2

            %velocity
            current_distance = pdist([p1; p2])/(t2-t1);
            dists(i+1) = current_distance;

            %if the current velocity is too big
            if current_distance > too_big

                %note that point (and the next 4) should be deleted (index for later)
                if length(eptrials(:,2))-2-i > guilty_by_association
                    deletions(i:i+guilty_by_association) = 1;
                end

                %move to the next point, but keep the first of the adjacent pair
                p2 = eptrials(i+2, 2:3);
                t2 = eptrials(i+2, 1);

                %each time it's too big, increase what is considered "too big"
                count = count + cum_mod;
                too_big = too_big + count;

            %if it's not too big
            else

                %reset what is considered "too big"
                too_big = too_big_origin;
                count = 0;

                %update points
                p1 = eptrials(i+1, 2:3);
                p2 = eptrials(i+2, 2:3);
                t1 = eptrials(i+1, 1);
                t2 = eptrials(i+2, 1);

            end
        end

        %index to delete dubious points
        deletions = logical(deletions);
        eptrials(deletions, 2:3) = NaN;

        %replace deleted points with interpolated values
        non_nan_pos = ~isnan(eptrials(:,2)) & ~isnan(eptrials(:,3)); %index
        eptrials(:,2) = interp1(eptrials(non_nan_pos, 1), eptrials(non_nan_pos, 2), eptrials(:,1), 'linear');
        eptrials(:,3) = interp1(eptrials(non_nan_pos, 1), eptrials(non_nan_pos, 3), eptrials(:,1), 'linear');

        %resample time and interpolate new position values
        pos(:,1) = (0:0.01:max(eptrials(:,1)))';
        pos(:,2) = interp1(eptrials(:,1), eptrials(:,2), pos(:,1), 'linear');
        pos(:,3) = interp1(eptrials(:,1), eptrials(:,3), pos(:,1), 'linear');
        eptrials = pos;
        
        %figure; plot(eptrials(:,2), eptrials(:,3))
    end

%add spike events (column 4)
eptrials = spk_evts(eptrials, ttun, spiketimes, initial_clock);
    function eptrials = spk_evts(eptrials, ttun, spiketimes, initial_clock)
    % spk_evts adds a column of spike events to eptrials. 
    %
    %   Spike events are coded as tetrode_number.cluster_number (e.g., 8.02 
    %       indicates the 2nd cluster on the 8th tetrode.
    %   Video re(samples) are coded as a 1. This helps with determining
    %       time by allowing you to sum the number of (equally spaced) time
    %       samples the occur during some period.
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   ttun= the tetrode and unit number in rows
    %   spiketimes = the spike timestamps for the tetrode and unit from ttun
    %   initial_clock = neuralynx start time
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %

        %add column
        eptrials = [eptrials ones(length(eptrials),1)];

        %some sessions had no clusters...
        if ~isempty(ttun)
            %build event_times, a two column matrix of times and cluster_ids
            event_times = [];
            for cluster = 1:size(ttun,1)
                cluster_id = ttun(cluster,1) + ttun(cluster,2)/100; %(tetrode.cluster)
                st_temp = [correct_time(spiketimes{cluster}, initial_clock)' repmat(cluster_id, length(spiketimes{cluster}),1)];
                event_times = [event_times; st_temp];
            end

            %combine eptrials and event_times
            event_times = [event_times(:,1) nan(length(event_times),2) event_times(:,2)];
            eptrials = sortrows([eptrials; event_times]);
        end
    end

%add flag events (column 5)
eptrials = flag_evts(eptrials, REWARD_L, REWARD_R, initial_clock);
    function eptrials = flag_evts(eptrials, REWARD_L, REWARD_R, initial_clock)
    % flag_evts adds a column of flag events to eptrials. flag_evts were
    % manually coded by Mark, an undergraduate student of Lindsey's. They
    % are double checked for errors (internal consistency) by a later function.
    %
    %   Flag events are coded as follows...
    %       Left reward instant = 1
    %       Right reward instant = 2
    %
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   REWARD_L = Video flag of left rewards
    %   REWARD_R = Video flag of right rewards
    %   initial_clock = the uncorrected session start time
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
    
        %add column
        eptrials = [eptrials zeros(length(eptrials),1)];

        %build matrices from flags
        RL_flags = [correct_time(REWARD_L, initial_clock)' nan(length(REWARD_L),3) ones(length(REWARD_L),1)];
        RR_flags = [correct_time(REWARD_R, initial_clock)' nan(length(REWARD_R),3) repmat(2, size(REWARD_R'))];

        %combine eptrials and flags
        eptrials = sortrows([eptrials; RL_flags; RR_flags]);
    end    

%delete events before the first and after the last position points
eptrials(eptrials(:,1) < min(eptrials(~isnan(eptrials(:,2)),1)), :) = [];
eptrials(eptrials(:,1) > max(eptrials(~isnan(eptrials(:,2)),1)), :) = [];

%interpolate position data at event times
eptrials = interp_at_nans(eptrials);

%add trial numbers (column 6)
eptrials = trial_nums(eptrials, TRIAL_START, initial_clock);
    function eptrials = trial_nums(eptrials, TRIAL_START, initial_clock)
    % trial_nums adds a column of trial numbers to eptrials. 
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   TRIAL_START = start times of each trial
    %   initial_clock = the uncorrected session start time
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
    
        %add column
        eptrials = [eptrials nan(length(eptrials),1)];

        %import times

        trial_starts = [TRIAL_START' nan(length(TRIAL_START),4) ones(length(TRIAL_START),1)];

        %correct times
        trial_starts(:,1) = (trial_starts(:,1) - repmat(initial_clock, size(trial_starts(:,1))))./1000000;

        %sort eptrials and trial_start flags
        eptrials = sortrows([eptrials; trial_starts]);

        %find begining and end of every trial in terms of rows
        trialstart_rows = find(eptrials(:,6)==1);
        trialend_rows = [trialstart_rows - ones(size(trialstart_rows)); length(eptrials)];
        trialend_rows(1) = [];

        %set row numbers
        for trial = 1:length(TRIAL_START)
            eptrials(trialstart_rows(trial):trialend_rows(trial),6) = trial;
        end
    end

%add choice type (column 8)
eptrials = LR_choice_types(eptrials);
    function eptrials = LR_choice_types(eptrials)
    % LR_trial_types adds a column of trial type to eptrials. This
    % indicates which light was on (i.e., where the reward was located).
    %
    %   Trial types are coded as follows...
    %       Left light was on = 1
    %       Right light was on = 2
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   trial_numbers = a vector of trial numbers (typically 1:40)
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %
    
        %add column
        eptrials = [eptrials nan(length(eptrials),2)];

        %left and right trial numbers
        left_trials = unique(eptrials(eptrials(:,5)==1,6));
        right_trials = unique(eptrials(eptrials(:,5)==2,6));

        %set rows
        eptrials(ismember(eptrials(:,6), left_trials), 8) = 1;
        eptrials(ismember(eptrials(:,6), right_trials), 8) = 2;
    end

%add trial accuracy (column 9)
eptrials = trial_accuracy(eptrials);
    function eptrials = trial_accuracy(eptrials)
    % trial_accuracy adds a column of choice accuracy type to eptrials. 
    % This indicates whether the rat turned in the same direction as the 
    % light (ie. correct or error).
    %
    %   Trial accuracy is coded as follows...
    %       Correct = 1
    %
    
        %add column
        eptrials = [eptrials ones(length(eptrials),1)];

    end

%verify reward/error flags
eptrials = verify_rwd_flg(eptrials);
    function eptrials = verify_rwd_flg(eptrials)
    %retrospectively verifies accuracy of reward and error flags, and then
    %corrects the inaccurate flags and reports them in command window. It
    %does this by making sure that the location of the rat at each reward
    %instant is closest to the mean location of the associated reward type.
    %

        %find x distance between location of rat at flag instant and mean 
        %of flag instances
        L_dists_to_L_mean = eptrials(ismember(eptrials(:,5),[1 3]), 2) - repmat(mean(eptrials(ismember(eptrials(:,5),[1 3]), 2)), length(eptrials(ismember(eptrials(:,5),[1 3]), 2)),1);
        L_dists_to_R_mean = eptrials(ismember(eptrials(:,5),[1 3]), 2) - repmat(mean(eptrials(ismember(eptrials(:,5),[2 4]), 2)), length(eptrials(ismember(eptrials(:,5),[1 3]), 2)),1);
        R_dists_to_L_mean = eptrials(ismember(eptrials(:,5),[2 4]), 2) - repmat(mean(eptrials(ismember(eptrials(:,5),[1 3]), 2)), length(eptrials(ismember(eptrials(:,5),[2 4]), 2)),1);
        R_dists_to_R_mean = eptrials(ismember(eptrials(:,5),[2 4]), 2) - repmat(mean(eptrials(ismember(eptrials(:,5),[2 4]), 2)), length(eptrials(ismember(eptrials(:,5),[2 4]), 2)),1);

        %flagging error if rat was closer to the other side
        L_incorrect_flag = (abs(L_dists_to_R_mean) - abs(L_dists_to_L_mean))<0;
        R_incorrect_flag = (abs(R_dists_to_L_mean) - abs(R_dists_to_R_mean))<0;

        %times associated with flags (and flags)
        l_rwd_times = eptrials(ismember(eptrials(:,5),[1 3]),[1 5]);
        r_rwd_times = eptrials(ismember(eptrials(:,5),[2 4]),[1 5]);

        %count errors
        rwd_flag_errors = 0;

        %correct any errors
        if sum(double(L_incorrect_flag))>0
            %count
            rwd_flag_errors = rwd_flag_errors + sum(double(L_incorrect_flag))>0;
            
            if ~isempty(l_rwd_times(L_incorrect_flag & l_rwd_times(:,2)==1,1))
                eptrials(ismember(eptrials(:,1), l_rwd_times(L_incorrect_flag & l_rwd_times(:,2)==1, 1)) & eptrials(:,5) > 0, 5) = 2;
            end
            if ~isempty(l_rwd_times(L_incorrect_flag & l_rwd_times(:,2)==3,1))
                eptrials(ismember(eptrials(:,1), l_rwd_times(L_incorrect_flag & l_rwd_times(:,2)==3, 1)) & eptrials(:,5) > 0, 5) = 4;
            end
        end
        if sum(double(R_incorrect_flag))>0
            %count
            rwd_flag_errors = rwd_flag_errors + sum(double(R_incorrect_flag))>0;
            
            if ~isempty(r_rwd_times(R_incorrect_flag & r_rwd_times(:,2)==2,1))
                eptrials(ismember(eptrials(:,1), r_rwd_times(R_incorrect_flag & r_rwd_times(:,2)==2, 1))  & eptrials(:,5) > 0, 5) = 1;
            end
            if ~isempty(r_rwd_times(R_incorrect_flag & r_rwd_times(:,2)==4,1))
                eptrials(ismember(eptrials(:,1), r_rwd_times(R_incorrect_flag & r_rwd_times(:,2)==4, 1))  & eptrials(:,5) > 0, 5) = 3;
            end
        end

        %report number of errors corrected
        if rwd_flag_errors > 0
            display(strcat([num2str(rwd_flag_errors), ' misflagged rewards/errors were corrected']))
        end
    end

%establish distance units (cm) by comparing distance between reward
%locations with distance measured on the actual maze
avgXYrwdL = [mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==1,3))];
avgXYrwdR = [mean(eptrials(eptrials(:,5)==2,2)) mean(eptrials(eptrials(:,5)==2,3))];
rwd_dist = pdist([avgXYrwdL; avgXYrwdR]);
cm_per_matlab_unit = 115.5/rwd_dist;

%coordinates of maze center
comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;

%correct maze orientation
eptrials = rotate_cntr(eptrials, comx, comy, avgXYrwdL, avgXYrwdR);
    function eptrials = rotate_cntr(eptrials, comx, comy, avgXYrwdL, avgXYrwdR)
    % rotate_cntr attempts to correct for angular misplacement of the maze 
    % below the camera by rotating every rat position around comxy with the
    % goal of leveling out the reward locations
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   comx = x coordinate of maze center
    %   comy = y coordinate of maze center
    %   avgXYrwdR = the average x position and y position [xpos
    %       ypos] of the rat at all right reward flags
    %   avgXYrwdL = the average x position and y position [xpos
    %       ypos] of the rat at all left reward flags
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    %

        %calculate the rotation angle required to level out the rewards
        rotang = rotation_angle(avgXYrwdL, avgXYrwdR, [comx comy]);
            function rotang = rotation_angle(avgXYrwdR, avgXYrwdL, comxy)
            %rotation_angle calculates rotang, the rotation angle required to 
            %level out the reward locations. Reward location inputs are
            %reversed to match lindsey/mark's manual coding.
            %
            % INPUTS
            %   avgXYrwdR = the average x position and y position [xpos
            %       ypos] of the rat at all right reward flags
            %   avgXYrwdL = the average x position and y position [xpos
            %       ypos] of the rat at all left reward flags
            %   comxy = [x_loc_centerofmaze y_loc_centerofmaze]
            %
            % OUTPUTS
            %   rotang = the rotation angle in clockwise degrees
            %

                %calculating what the new right reward location should
                    %set the reward location to the origin before finding the angle
                    oldX_R_origin = avgXYrwdR(1) - comxy(1);
                    oldY_R_origin = avgXYrwdR(2) - comxy(2);

                    %radius can be found using pythag theorum
                    radius = sqrt(oldX_R_origin^2 + oldY_R_origin^2);

                    %new y rwd coordinate is the mean y coordinate of both rwds
                    newY_R_origin = oldY_R_origin + (avgXYrwdL(2) - avgXYrwdR(2))/2;

                    %the new x coordinate can be found using pythag theorum
                    newX_R_origin = sqrt(radius^2 - newY_R_origin^2);

                %determine rotation angle (in clockwise degrees)
                    %clockwise turn if the new R_RWD y is lower than the old
                    if newY_R_origin < oldY_R_origin
                        rotang = acosd((2*radius^2 - pdist([oldX_R_origin oldY_R_origin; newX_R_origin newY_R_origin], 'euclidean')^2)/(2*radius^2));

                    %counterclickwise if the new R_RWD y coordinate is higher
                    else
                        rotang = -acosd((2*radius^2 - pdist([oldX_R_origin oldY_R_origin; newX_R_origin newY_R_origin], 'euclidean')^2)/(2*radius^2));
                    end
            end
            
        %verify rotation angle
        if abs(rotang) > 20
            error('rotation angle implausible. check rotation angle')
        end

        %report rotation angle (negative for lindsey)
        display(strcat(['X,Y coordinates were rotated clockwise ', num2str(-rotang), ' degrees']))

        %rotate all x y coordinates
        newpts = rotate_pts(rotang*-1, [eptrials(:,2) eptrials(:,3)], [comx comy]);
            function newpts = rotate_pts(rotang, pts, center)
            %rotate_pts takes a list of i j coordinates "pts" and returns them 
            %rotated around "center" by the angle "rotang"
            %
            % INPUTS
            %   rotang = rotation angle counter clockwise in degrees
            %   pts = points to rotate
            %   center = point you would like to rotate the pts about 
            %
            % OUTPUTS
            %   newpts = list of rotated i j coordinates size(2,length(pts))
            %

                %translate to radians
                rotang = rotang/57.2957795;

                %build rotation matrix (some calc II thing)
                romat=[cos(rotang) -sin(rotang);sin(rotang) cos(rotang)];

                %rotation occurs around origin, so we temporarily center all points around the origin
                pts(:,1)=pts(:,1)-center(:,1);
                pts(:,2)=pts(:,2)-center(:,2);

                %apply rotation
                newpts=(romat*pts')';

                %undo centering
                newpts(:,1)=newpts(:,1)+center(:,1);
                newpts(:,2)=newpts(:,2)+center(:,2);

            end

        %rewrite the x and y coordinates of eptrials
        eptrials(:,2) = newpts(:,1);
        eptrials(:,3) = newpts(:,2);

    end

%set center of maze to (1000,1000)
comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;
x_correction = 1000-comx;
y_correction = 1000-comy;
eptrials(:,2) = eptrials(:,2) + repmat(x_correction, size(eptrials(:,2)));
eptrials(:,3) = eptrials(:,3) + repmat(y_correction, size(eptrials(:,3)));
comx = mean([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==2,2))]);
comy = mean([mean(eptrials(eptrials(:,5)==1,3)) mean(eptrials(eptrials(:,5)==2,3))])-100;

%establish maze sections (column 10)
[eptrials, strt, strt_idx, not_start, section_boundaries] = maze_sections(eptrials, comx, comy);
    function [eptrials, strt, strt_idx, not_start, section_boundaries] = maze_sections(eptrials, comx, comy)
            %maze_sections defines and builds the maze sections around the
            %center of the maze as defined by comx/y. It does not deal with the
            %iti section.
            %
            % INPUTS
            %   eptrials = the overall output; a matrix of all relevant information
            %   comx = x coordinate of maze center
            %   comy = y coordinate of maze center
            %
            % OUTPUTS
            %   eptrials = the overall output; a matrix of all relevant information
            %   strt = maze section boundaries of the start section
            %   strt_idx = index for x y points that fall within strt   
            %   not_strt = ~strt_idx
            %

                %add column
                eptrials = [eptrials nan(length(eptrials),1)];

                %establishes maze section boundaries [xlow xhigh ylow yhigh]
                strt = [comx-50 comx+50  comy-200 comy-90]; %start area
                stem = [comx-30 comx+30 comy-90 comy+60]; %common stem
                chce = [comx-45 comx+45 comy+60 comy+155]; %choice area
                chmL = [comx+45 comx+135 comy+75 comy+155]; %approach arm left
                chmR = [comx-135 comx-45 comy+75 comy+155]; %approach arm right
                rwdL = [comx+135 comx+200 comy+60 comy+155]; %reward area left
                rwdR = [comx-200 comx-135 comy+60 comy+155]; %reward area right
                
                section_boundaries = [strt;stem;chce;chmL;chmR;rwdL;rwdR];
                
                %anti-start index
                not_start = eptrials(:,2)<strt(1) | eptrials(:,2)>strt(2) | eptrials(:,3)<strt(3) | eptrials(:,3)>strt(4);

                %section indices
                strt_idx = eptrials(:,2)>=strt(1) & eptrials(:,2)<=strt(2) & eptrials(:,3)>=strt(3) & eptrials(:,3)<=strt(4);
                stem_idx = eptrials(:,2)>=stem(1) & eptrials(:,2)<=stem(2) & eptrials(:,3)>=stem(3) & eptrials(:,3)<=stem(4);
                chce_idx = eptrials(:,2)>=chce(1) & eptrials(:,2)<=chce(2) & eptrials(:,3)>=chce(3) & eptrials(:,3)<=chce(4);
                chmL_idx = eptrials(:,2)>=chmL(1) & eptrials(:,2)<=chmL(2) & eptrials(:,3)>=chmL(3) & eptrials(:,3)<=chmL(4);
                chmR_idx = eptrials(:,2)>=chmR(1) & eptrials(:,2)<=chmR(2) & eptrials(:,3)>=chmR(3) & eptrials(:,3)<=chmR(4);
                rwdL_idx = eptrials(:,2)>=rwdL(1) & eptrials(:,2)<=rwdL(2) & eptrials(:,3)>=rwdL(3) & eptrials(:,3)<=rwdL(4);
                rwdR_idx = eptrials(:,2)>=rwdR(1) & eptrials(:,2)<=rwdR(2) & eptrials(:,3)>=rwdR(3) & eptrials(:,3)<=rwdR(4);

                %set rows for maze
                eptrials(strt_idx, 10) = 1;
                eptrials(stem_idx, 10) = 2;
                eptrials(chce_idx, 10) = 3;
                eptrials(chmL_idx, 10) = 4;
                eptrials(chmR_idx, 10) = 5;
                eptrials(rwdL_idx, 10) = 6;
                eptrials(rwdR_idx, 10) = 7;
    end

%add iti tag (column 10)
[eptrials iti] = set_iti(eptrials, comx, comy, ITI_ON, ITI_OFF, initial_clock, strt_idx, not_start);
    function [eptrials, iti, iti_idx, not_iti] = set_iti(eptrials, comx, comy, ITI_ON, ITI_OFF, initial_clock, strt_idx, not_start)
    % set_iti overwrites elements in column 10 (maze section) of eptrials 
    % to indicate that the intertrial interval is currently active. This 
    % is related to but morespecific than information about whether the 
    % rat is in the iti section relative to the center of the maze. It 
    % relies on manually flagging of the iti begining and end times. 
    % Because some of these manual flags are incorrect, set_iti verifies 
    % them, corrects them (or at least improves upon them) and reports 
    % errors in the command window
    %
    % INPUTS
    %   eptrials = the overall output; a matrix of all relevant information
    %   comx = x coordinate of maze center
    %   comy = y coordinate of maze center
    %   ITI_ON = first timestamp of intertrial interval
    %   ITI_OFF = last timestamp of intertrial interval
    %   initial_clock = the uncorrected session start time
    %   strt_idx = index for x y points that fall within strt   
    %   not_strt = ~strt_idx
    %
    % OUTPUTS
    %   eptrials = the overall output; a matrix of all relevant information
    %   not_iti = index for x y points outside of iti maze section
    %              


    	%iti section ambit and indices
     	iti = [comx+30 comx+180  comy-170 comy-40]; %iti section
     	iti_idx = eptrials(:,2)>=iti(1) & eptrials(:,2)<=iti(2) & eptrials(:,3)>=iti(3) & eptrials(:,3)<=iti(4);
     	not_iti = eptrials(:,2)<iti(1) | eptrials(:,2)>iti(2) | eptrials(:,3)<iti(3) | eptrials(:,3)>iti(4);

       	%set rows for iti
     	iti_starts = correct_time(ITI_ON, initial_clock);
      	iti_ends = correct_time(ITI_OFF, initial_clock);
      	iti_starts_ends = [iti_starts' iti_ends'];
       	iti_flag_error_count = 0;
      	iti_flag_error_trials = [];

      	%for all trials
      	for trial = 1:max(eptrials(:,6))

          	%find range by finding rows with closest times
           	[~,start_row] = min(abs(eptrials(:,1) - repmat(iti_starts_ends(trial, 1), size(eptrials(:,1)))));
           	[~,end_row] = min(abs(eptrials(:,1) - repmat(iti_starts_ends(trial, 2), size(eptrials(:,1)))));

           	%set rows in range
          	eptrials(start_row:end_row, 10) = 8;

           	%check for iti flagging error
               	%if during iti, rat enters start area outside of iti 
              	%bounds and stays there for 0.5s, then assume iti has 
               	%ended

               	%index from range
              	ep_range = zeros(length(eptrials(:,10)), 1);
              	ep_range(start_row:end_row) = 1;
              	ep_range = logical(ep_range);

              	%find vid_sample points during iti range that are located in start 
              	%area but not located in iti
               	oob_iti = ep_range & strt_idx & not_iti;
               	vid_oob_iti = oob_iti(eptrials(:,4)==1);

              	%preallocate iti exit time index
               	vid_iti_exits = zeros(size(vid_oob_iti));

               	%find first point of every 5 consecutive 1's
              	for win = 1:length(vid_oob_iti)-4
                  	if sum(vid_oob_iti(win:win+4)) == 5
                      	vid_iti_exits(win) = 1;
                       	iti_flag_error_count = iti_flag_error_count + 1;
                       	iti_flag_error_trials = [iti_flag_error_trials trial];
                       	break
                    end
                end

              	%find row of certain iti exit
               	iti_exits = zeros(size(oob_iti));
              	iti_exits(eptrials(:,4)==1) = vid_iti_exits; %resize
               	exit_row = find(iti_exits==1, 1);

             	%use slightly earlier iti exit (last time point in iti area)
              	if ~isempty(exit_row)
                 	idx_last_iti = max(eptrials(iti_idx & not_start & eptrials(:,1)<eptrials(exit_row,1), 1));
                 	exit_row = find(eptrials(:,1) == idx_last_iti, 1);
                end

                 	%correct eptrials rows
                  	eptrials(exit_row:end_row, 10) = NaN;
                    
        end

      	%report iti flag errors
      	if iti_flag_error_count > 0
            display(strcat([num2str(iti_flag_error_count), ' misflagged iti periods were improved']))
          	%iti_flag_error_trials
      	end
    end

%redfine iti section area based on an average of the iti positions
[~, iti_idx, not_iti, section_boundaries] = shift_iti(eptrials, iti, section_boundaries);
    function [iti, iti_idx, not_iti, section_boundaries] = shift_iti(eptrials, iti, section_boundaries)
    %locates the iti maze section. Good locating of this
    %section enables the deletion of bad points, and the detection of bad
    %iti flags.

        %preallocate
        max_x = nan(max(eptrials(:,6)),1);
        min_x = nan(max(eptrials(:,6)),1);
        max_y = nan(max(eptrials(:,6)),1);
        min_y = nan(max(eptrials(:,6)),1);
        
        %iterate to find trial variation
        for trial = 1:max(eptrials(:,6))
            try
                max_x(trial) = max(eptrials(eptrials(:,6)==trial & eptrials(:,10)==8, 2));
                min_x(trial) = min(eptrials(eptrials(:,6)==trial & eptrials(:,10)==8, 2));
                max_y(trial) = max(eptrials(eptrials(:,6)==trial & eptrials(:,10)==8, 3));
                min_y(trial) = min(eptrials(eptrials(:,6)==trial & eptrials(:,10)==8, 3));
            catch
            end
        end
        
        %take medians
        minmax_xy = nanmedian([min_x max_x min_y max_y]);
        
        %center iti
        mean_pos = [mean([minmax_xy(1) minmax_xy(2)]) mean([minmax_xy(3) minmax_xy(4)])];
        
        %current iti
        mean_current_iti = [mean(iti(1:2)) mean(iti(3:4))];
        
        %difference
        pos_dif = mean_pos - mean_current_iti;
        
        %update iti position
        iti(1:2) = iti(1:2) + [pos_dif(1) pos_dif(1)];
        iti(3:4) = iti(3:4) + [pos_dif(2) pos_dif(2)];
        
        %update iti indices
        iti_idx = eptrials(:,2)>=iti(1) & eptrials(:,2)<=iti(2) & eptrials(:,3)>=iti(3) & eptrials(:,3)<=iti(4);
     	not_iti = eptrials(:,2)<iti(1) | eptrials(:,2)>iti(2) | eptrials(:,3)<iti(3) | eptrials(:,3)>iti(4);
        
        %update sections
        section_boundaries = [section_boundaries;iti];

    end

%add carry tag (column 10)
eptrials = set_carry(eptrials, cm_per_matlab_unit, strt, not_start, iti_idx);
    function eptrials = set_carry(eptrials, cm_per_matlab_unit, strt, not_start, iti_idx)
        %set_carry overwrites elements in column 10 (maze section) of 
        %eptrials to indicate that the rat is being carried. It is based on
        %sudden changes in velocity just after the events that typically
        %proceed pick-ups: reward and the end of the ITI.
        %
        % INPUTS
        %   eptrials = the overall output; a matrix of all relevant information
        %   cm_per_matlab_unit = conversion between xy coordniate space and real-life distance
        %   strt = maze section boundaries of the start section
        %   not_strt = index of coordinate position outside of start
        %
        % OUTPUT
        %   eptrials = the overall output; a matrix of all relevant information
        % 
        

        %set rows for carry lable (i.e., 'not maze'; will appropriately overwrite some of the above maze categorizations)
        for trial = 1:max(eptrials(:,6))

            %one trial of eptrials
            local_eptrials = eptrials(eptrials(:,6)==trial, :);
            loc_strt_idx = local_eptrials(:,2)>=strt(1) & local_eptrials(:,2)<=strt(2) & local_eptrials(:,3)>=strt(3) & local_eptrials(:,3)<=strt(4);

            %calulate smoothed velocity
            	%just vid samples
                vid_local_eptrials = local_eptrials(local_eptrials(:,4)==1, 1:3);

              	%preallocate
              	veloc = zeros(length(vid_local_eptrials(:,1)-1),1);

              	%place and time
              	p1 = vid_local_eptrials(1:end-1, 2:3);
             	p2 = vid_local_eptrials(2:end, 2:3);
              	t1 = vid_local_eptrials(1:end-1, 1);
               	t2 = vid_local_eptrials(2:end, 1);

              	%calculate
               	for instant = 1:length(vid_local_eptrials(:,1))-1
                 	veloc(instant+1) = pdist([p1(instant,:); p2(instant,:)])/(t2(instant)-t1(instant));
                end

               	veloc = veloc.*cm_per_matlab_unit;%correct to distance/seconds
               	veloc = smooth(veloc, 50);%smooth
               	local_velocity = [vid_local_eptrials(:,1) veloc];%output
               	high_veloc = 35;
               	low_veloc = 25;

            %pick up after reward
                rwd_time_allotment = 1; %arbitrary
                post_rwd = local_eptrials(local_eptrials(:,5)>0,1) + rwd_time_allotment; %earliest possible end of reward
                hiveloc_instants = local_velocity(local_velocity(:,2)>high_veloc,1);%high velocity
                rwd_pick_up = min(hiveloc_instants(hiveloc_instants > post_rwd)); %begin carry 
                [~,start_row] = min(abs(eptrials(:,1) - repmat(rwd_pick_up, size(eptrials(:,1))))); %find row

            %put down in start area
                last_iti = max(local_eptrials(local_eptrials(:,10)==8,1)); %end of iti
                
                %if session ends during iti
                if isempty(last_iti)
                    last_iti = max(local_eptrials(:,1));
                end
                
                loveloc_instants = local_velocity(local_velocity(:,2)<low_veloc,1);%low velocity

                %if rat enters start area after leaving iti and decelerates
                if ~isempty(loveloc_instants(loveloc_instants>last_iti)) && ~isempty(local_eptrials(local_eptrials(:,10)==1 & local_eptrials(:,1)>last_iti, 1))

                    start_put_down = min(loveloc_instants(loveloc_instants > last_iti)); %end carry

                    %make sure put down is after the first start area after the maximum iti
                    if start_put_down > min(local_eptrials(loc_strt_idx & local_eptrials(:,1)>last_iti,1))

                        [~,end_row] = min(abs(eptrials(:,1) - repmat(start_put_down, size(eptrials(:,1))))); %find row

                    % if it isnt, then default to trial end (velocity trick failed) 
                    else
                        end_row = find(eptrials(:,1)==local_eptrials(end,1),1);%go to end of session
                    end

                %else, use end of trial
                else
                    end_row = find(eptrials(:,1)==local_eptrials(end,1),1);%go to end of session
                end
                
            %set rows
            ep_range = zeros(length(eptrials(:,10)), 1);
            ep_range(start_row:end_row) = 1;
            ep_range = logical(ep_range);
            eptrials(ep_range & eptrials(:,10)~=8, 10) = 9;
            
            %santity check - does iti end in iti?  
            %
            %time of last iti tag
            tot = max(eptrials(eptrials(:,6)==trial & eptrials(:,10) == 8, 1));
            %last in iti section
            toi = max(eptrials(eptrials(:,6)==trial & iti_idx, 1));
            %fix
            if tot > toi
                eptrials(eptrials(:,1)>=toi & eptrials(:,1)<=tot, 10) = 9;
            end

            %delete points after start_put_down that are out of start area
            eptrials(eptrials(:,6)==trial & eptrials(:,1)>eptrials(end_row,1) & not_start, 2:3) = NaN;     
        end
    end

%scrub location data in terms of above-defined maze sections
eptrials(isnan(eptrials(:,10)), 2:3) = NaN; %positions without sections
eptrials(eptrials(:,10) == 8 & not_iti, 2:3) = NaN;%miscategorized ITI

%delete events before the first and after the last position points
eptrials(eptrials(:,1) < min(eptrials(~isnan(eptrials(:,2)),1)), :) = [];
eptrials(eptrials(:,1) > max(eptrials(~isnan(eptrials(:,2)),1)), :) = [];

%interpolate position at scrubbed location data points
    if ~isempty(eptrials(isnan(eptrials(:,2)) & isnan(eptrials(:,3)),1))
        eptrials = interp_at_nans(eptrials);
    end

%smooth position data
eptrials(:,2) = smooth(eptrials(:,2), 10);
eptrials(:,3) = smooth(eptrials(:,3), 10);

%establish stem sections (column 11)
[eptrials, lost_trials] = stem_sections(eptrials, comx, comy, lost_trials);
    function [eptrials, lost_trials] = stem_sections(eptrials, comx, comy, lost_trials)
            %stem_sections defines and builds the stem sections around the
            %center of the maze as defined by comx/y. They are intended to fit 
            %neatly within the bounds of the stem as defined by stem sections.
            %Importantly, the stem section label is only applied to timepoints
            %after and before the rat's true entrance to and exit from the
            %stem on that trial, as defined below.
            %
            % INPUTS
            %   eptrials = the overall output; a matrix of all relevant information
            %   comx = x coordinate of maze center
            %   comy = y coordinate of maze center
            %
            % OUTPUTS
            %   eptrials = the overall output; a matrix of all relevant information
            %       

            %add column
            eptrials = [eptrials nan(length(eptrials),1)];

            %stem sections
            stem1 = [comx-30 comx+30 comy-90 comy-52.5];
            stem2 = [comx-30 comx+30 comy-52.5 comy-15];
            stem3 = [comx-30 comx+30 comy-15 comy+22.5];
            stem4 = [comx-30 comx+30 comy+22.5 comy+60];

            %set rows for sub-stem sections
            for trial = 1:max(eptrials(:,6))

                %the last time point in the start area on this trial before reward
                enter_stem = max(eptrials(eptrials(:,6)==trial & eptrials(:,10)==1 & eptrials(:,1) < min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,5), 1:4), 1)) ,1));
                %the first time point in the choice area (or beyond) on this trial 
                exit_stem = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,10), [3 4 5 6 7]) & eptrials(:,1) > enter_stem, 1));
                
                if ~isempty(enter_stem) && ~isempty(exit_stem)

                    eptrials(eptrials(:,1)>enter_stem & eptrials(:,1)<exit_stem & eptrials(:,3)>stem1(3) & eptrials(:,3)<stem1(4),11) = 1;
                    eptrials(eptrials(:,1)>enter_stem & eptrials(:,1)<exit_stem & eptrials(:,3)>stem2(3) & eptrials(:,3)<stem2(4),11) = 2;
                    eptrials(eptrials(:,1)>enter_stem & eptrials(:,1)<exit_stem & eptrials(:,3)>stem3(3) & eptrials(:,3)<stem3(4),11) = 3;
                    eptrials(eptrials(:,1)>enter_stem & eptrials(:,1)<exit_stem & eptrials(:,3)>stem4(3) & eptrials(:,3)<stem4(4),11) = 4;
                    
                else
                    lost_trials = [lost_trials; trial];
                    warning('Unable to define stem entrance and exit')
                end

            end
    end

%add velocity (column 12)
eptrials = vid_velocity(eptrials, cm_per_matlab_unit);
    function eptrials = vid_velocity(eptrials, cm_per_matlab_unit)
    %vid_velocity adds a column with the estimated instantaneous
    %velocity at every video sample time-point (row). Other rows are nan.
    %
    % INPUTS
    %   eptrials = the overall output; a matrix of all relevant information
    %   cm_per_matlab_unit = conversion between xy coordniate space and real-life distance
    %
    % OUTPUT
    %   eptrials = the overall output; a matrix of all relevant information
    %    
        
        %add column
        eptrials = [eptrials nan(length(eptrials),1)];

        %calculate distances between every video sample
        eptrials_vid = eptrials(eptrials(:,4) == 1, 2:3);

        pos1 = [eptrials_vid(1:end-1,1)'; eptrials_vid(1:end-1,2)'];
        pos2 = [eptrials_vid(2:end,1)'; eptrials_vid(2:end,2)'];

        distances = zeros(length(eptrials_vid),1);
        for i = 1:length(eptrials_vid)-1
            distance = dist([pos1(:,i) pos2(:,i)]);
            distances(i+1) = distance(1,2);
        end

        %calulate velocity from dists_overall (see above)
        velocity = (distances.*cm_per_matlab_unit)./0.01;

        %set rows
        eptrials(eptrials(:,4) == 1,12) = smooth(velocity, 40);
        
    end

%add light on/off (column 13)
eptrials = light_status(eptrials, initial_clock, LIGHT_ON, FAKE_ON, LIGHT_OFF, TRIAL_START);
    function eptrials = light_status(eptrials, initial_clock, LIGHT_ON, FAKE_ON, LIGHT_OFF, TRIAL_START)
    % light_status adds a column to eptrials indicating whether a light is
    % on and which light, left or right
    %
    %   Light status is coded as follows...
    %       Light on = 1
    %       Fake light on = 3
    %
    % INPUT VARIABLES
    %   eptrials = the overall output; a matrix of all relevant information
    %   initial_clock = the uncorrected session start time
    %   LIGHT_ON = Timestamps for left light turning on
    %   LIGHT_OFF = Timestamps for left light turning off
    %
    % OUTPUT VARIABLE
    %   eptrials = the overall output; a matrix of all relevant information
    % 
        
        %add column
        eptrials = [eptrials nan(length(eptrials),1)];
            
            %set rows for light on
            Light_starts = (LIGHT_ON - repmat(initial_clock, size(LIGHT_ON)))./1000000; 
            Light_ends = (LIGHT_OFF - repmat(initial_clock, size(LIGHT_OFF)))./1000000;
            Light_starts_ends = [Light_starts' Light_ends'];    
            for Light = 1:size(Light_starts_ends,1)

                %find range
                [~,start_row] = min(abs(eptrials(:,1) - repmat(Light_starts_ends(Light, 1), size(eptrials(:,1)))));
                [~,end_row] = min(abs(eptrials(:,1) - repmat(Light_starts_ends(Light, 2), size(eptrials(:,1)))));

                %set rows in range
                eptrials(start_row:end_row, 13) = 1;
            end

            
            %set rows for fake light on
            FLight_starts = ((FAKE_ON - repmat(initial_clock, size(FAKE_ON)))./1000000)'; 
            %FLight_ends = (TRIAL_START - repmat(initial_clock, size(TRIAL_START)))./1000000;
            for Light = 1:size(FLight_starts,1)

                %find range
                [~,start_row] = min(abs(eptrials(:,1) - repmat(FLight_starts(Light, 1), size(eptrials(:,1)))));
                %[~,end_row] = min(abs(eptrials(:,1) - repmat(Light_starts_ends(Light, 2), size(eptrials(:,1)))));

                %set rows in range
                eptrials(start_row, 13) = 3;
            end
    end

%Remove extraneous points before and after session
    %find the first and last time points that recieved palce, trial, and
    %section categorizations
    start_trial_labels = max([min(eptrials(~isnan(eptrials(:,2)),1)) min(eptrials(~isnan(eptrials(:,6)),1)) min(eptrials(~isnan(eptrials(:,10)),1))]);
    end_trial_labels = min([max(eptrials(~isnan(eptrials(:,2)),1)) max(eptrials(~isnan(eptrials(:,6)),1)) max(eptrials(~isnan(eptrials(:,10)),1))]);

    %delete events before and after these time points
    eptrials(eptrials(:,1) < start_trial_labels, :) = [];
    eptrials(eptrials(:,1) > end_trial_labels, :) = [];

    %reset time to a 0s start
    eptrials(:,1) = eptrials(:,1) - repmat(min(eptrials(:,1)), size(eptrials(:,1)));
    
%REUSED SUBFUNCTIONS
function new_time = correct_time(old_time, initial_clock)
    %corrects neuralynx clock to a 0s session start time
	new_time = (old_time - repmat(initial_clock, size(old_time)))./1000000;
end
function eptrials = interp_at_nans(eptrials)
    %interpolate_at_nans finds missing position information in eptrials
    %and interpolates to fill
    %
  
    %index for missing position values
    real_x = eptrials(~isnan(eptrials(:,2)),2);
    real_y = eptrials(~isnan(eptrials(:,3)),3);
    all_time = eptrials(:,1);
    
    %only including unique time points among the missing rows
    [real_time_x, idx_x] = unique(eptrials(~isnan(eptrials(:,2)),1));
    [real_time_y, idx_y] = unique(eptrials(~isnan(eptrials(:,3)),1));
    
    %index real_x and real_y to match real_time unique elements
    eptrials(:,2) = interp1(real_time_x, real_x(idx_x), all_time);
    eptrials(:,3) = interp1(real_time_y, real_y(idx_y), all_time);  
    
end

end