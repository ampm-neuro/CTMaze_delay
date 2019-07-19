function [light_on_times, stem_ent_times, stem_exit_times, rwd_times] = ALL_eventtimes(stage_numbers)
%all performs some function copy pasted below on every session. typically
%for plotting


% stage numbers
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day

%set counters
%count = 0;


rate_matrix_4shuf = nan(42, 9, 360);

%min_pop = 8;
%popday = 0;

light_on_times = [];
stem_ent_times = [];
stem_exit_times = [];
rwd_times = [];

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);



sesh_num = 0;

%iterate through subjects...and stages...and sessions
for subject = 1:length_subjects
    
    
    %figure;hold on
    
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
        task = file_names_stages{:}(stage,1).name
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(file_names_stages{:}(stage,1).name, '/*.mat')));
        file_list_sessions(1:2) = [];
        length_sessions = size(file_list_sessions,1) -1;
    
        %iterate through sessions
        for session = 1:length_sessions
            
            
            if session <10
                session_num = strcat('0', num2str(session));
            else
                session_num = session;
            end
            
            %figure;hold on
            
            %load
            eptrials = [];
            origin_file = [];
            clusters = [];
            light_on_sect = [];
            out_reward = [];
            load(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/',num2str(rat),'/' ,num2str(task),'/' ,num2str(session_num), '.mat'));
            
            if isempty(clusters) %|| length(clusters(:,1)) < min_pop
                continue
            else
                %popday = popday+1;
            end    
               
            session

            
               %rejected
               all_trials = unique(eptrials(:,6));
               rejected_trials = reject_trials('all', lost_trials, out_stem, light_on_sect, [1 2.1 2.2 2.3 2.4], out_reward);
               included_trials = setdiff(all_trials,rejected_trials);
               error_trials = unique(eptrials(eptrials(:,9)==0,6));
               included_trials = setdiff(included_trials,error_trials);
               

               for trial = included_trials'
                                  
                   trial_start = min(eptrials(eptrials(:,6)==trial,1));
                   
         
                   %first light
                   fl = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13),[1 2]),1)) - trial_start;
                   if ~isempty(fl)
                        light_on_times = [light_on_times; fl];
                        
                        ent = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11),1),1)) - trial_start;
                        if ~isempty(ent)
                            stem_ent_times = [stem_ent_times; ent];
                        else
                            stem_ent_times = [stem_ent_times; nan];
                        end
                        
                        ext = max(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11),4),1)) - trial_start;
                        if ~isempty(ext)
                            stem_exit_times = [stem_exit_times; ext];
                        else
                            stem_exit_times = [stem_exit_times; nan];
                        end
                        
                        rwd = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,5),[1 2 3 4]),1)) - trial_start;
                        if ~isempty(rwd)
                            rwd_times = [rwd_times; rwd];
                        else
                            rwd_times = [rwd_times; nan];
                        end
                        
                   else
                       light_on_times = [light_on_times; nan];
                       stem_ent_times = [stem_ent_times; nan];
                       stem_exit_times = [stem_exit_times; nan];
                       rwd_times = [rwd_times; nan];
                        
                   end
                
               end
                %}
            
        end
    end
end
end
 