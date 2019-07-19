function ALL_traj_plot(stage_numbers, bw, aw, varargin)
%all performs some function copy pasted below on every session. typically
%for plotting

if nargin == 4
    prop_keep = varargin{1}
end

%set counters
figure
hold on

grn=[52 153 70]./255;
blu=[46 49 146]./255;
colors = [255 225 102; 255 178 102; 216.7500 82.8750 24.9900; 142 0 0; 102 205 255; 51 153 255;  0  113.9850  188.9550; 50 0 150];

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);

%iterate through subjects...and stages...and sessions
for subject = 1:length_subjects
    %print update
    rat = file_names_subjects{:}(subject,1).name
    %get all the things in subject folder...
    file_list_stages = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name)));
    file_list_stages(1:2) = [];
    file_names_stages = {file_list_stages([file_list_stages(:).isdir])};
    length_stages = size(file_names_stages{:},1);
    
    %iterate through stages...and sessions
    for stage = stage_numbers
        %print update
        task = file_names_stages{:}(stage,1).name;
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(file_names_stages{:}(stage,1).name, '/*.mat')));
        file_list_sessions(1:2) = [];
        length_sessions = size(file_list_sessions,1)-1;
    
        %iterate through sessions
        for session = 1:length_sessions
            
            if session <10
                session_num = strcat('0', num2str(session));
            else
                session_num = session;
            end
            
            eptrials = [];
            origin_file = [];
            clusters = [];
            light_on_sect = [];
            load(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/',num2str(rat),'/' ,num2str(task),'/' ,num2str(session_num), '.mat'));
                
            if session == 1
                sections(eptrials, section_boundaries)
            end
            
            for trial = 1:max(eptrials(:,6))
                if ~isempty(eptrials(eptrials(:,6)==trial & ~isnan(eptrials(:,13)), 1))
                    
                    
                    %moment = min(eptrials(eptrials(:,6)==trial, 1)); %start
                    moment = min(eptrials(eptrials(:,6)==trial & ~isnan(eptrials(:,13)), 1)); %light
                    %moment = max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4, 1)); %stem exit
                    %moment = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:, 5),1:2), 1)); %rwd
                    
                    if isempty(moment)
                        continue
                    end
                    
                    %moment
                    
                    %around; trialtype
                    %left_trl_idx = eptrials(:,6)==trial & eptrials(:,8)==1 & eptrials(:,1)>moment-1 & eptrials(:,1)<moment+1;
                    %right_trl_idx = eptrials(:,6)==trial & eptrials(:,8)==2 & eptrials(:,1)>moment-1 & eptrials(:,1)<moment+1;
                    
                    %before & after
                    
                    %before/after window
                    %bw = [-.5 0];
                    %aw = [0 2];
                    
                    if rand(1) < prop_keep
                    
                        before_trl_idx = eptrials(:,1)>moment+bw(1) & eptrials(:,1)<moment+bw(2);
                        after_trl_idx = eptrials(:,1)>moment+aw(1) & eptrials(:,1)<moment+aw(2);
                        
                        if mode(eptrials(before_trl_idx | after_trl_idx, 8)) == 1
                        
                            %plot(eptrials(before_trl_idx, 2), eptrials(before_trl_idx, 3), 'Color', 'k', 'LineWidth', 0.5, 'LineStyle', '-');
                            %plot(eptrials(after_trl_idx, 2), eptrials(after_trl_idx, 3), 'Color', colors(2,:)./255, 'LineWidth', 0.5, 'LineStyle', '-');
                            %plot(eptrials(eptrials(:,1)==moment, 2), eptrials(eptrials(:,1)==moment, 3), '.', 'Color', colors(2,:)./255, 'MarkerSize', 30);
                            
                            
                        elseif mode(eptrials(before_trl_idx | after_trl_idx, 8)) == 2
                            
                            plot(eptrials(after_trl_idx, 2), eptrials(after_trl_idx, 3), 'Color', colors(6,:)./255, 'LineWidth', 0.5, 'LineStyle', '-');
                            plot(eptrials(find(after_trl_idx==1, 1, 'first'), 2), eptrials(find(after_trl_idx==1, 1, 'first'), 3), '.', 'Color', colors(6,:)./255, 'MarkerSize', 30);

                            
                        end

                    end
                    
                end
            end
            
        end
    end
end
end






                
                
                
              