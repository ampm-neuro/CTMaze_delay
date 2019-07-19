function [comb_distmtx, flag_times] = ALL_lighton(stage_numbers)
%all performs some function copy pasted below on every session. typically
%for plotting


% stage numbers
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day

%set counters
%count = 0;


rate_matrix_4shuf = nan(42, 9, 360);

min_pop = 8;
popday = 0;

comb_distmtx = [];
flag_times = [];

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);



sesh_num = 0;

figure;hold on

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
        length_sessions = size(file_list_sessions,1);
    
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
            
            if isempty(clusters) || length(clusters(:,1)) < min_pop
                continue
            else
                popday = popday+1;
            end    
               
            session

            
               %rejected
               all_trials = unique(eptrials(:,6));
               rejected_trials = reject_trials('all', lost_trials, out_stem, light_on_sect, [1 2.1 2.2 2.3 2.4], out_reward);
               included_trials = setdiff(all_trials,rejected_trials);
               error_trials = unique(eptrials(eptrials(:,9)==0,6));
               included_trials = setdiff(included_trials,error_trials);
               
               
               %[~, stem_ent,stem_ext, rwd_time] = LR_dist_lite(eptrials, clusters, included_trials, .25);
               %comb_distmtx = [comb_distmtx; dist_mtx];
               
               %flag_times = [flag_times; [stem_ent stem_ext rwd_time]];
               
               
               %
               %plot maze sections
               sections(eptrials, section_boundaries)
               

               for trial = included_trials'
                                  
                   
         
                   %first light
                   fl = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13),[1 2]),1));
                   
                   %if ~isempty(fl)
                   
                       %first light x and y
                       flx = eptrials(eptrials(:,1) == fl, 2);
                       fly = eptrials(eptrials(:,1) == fl, 3);

                       %plot light on
                       if mode(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13),[1 2]), 13) == 1)
                           plot(flx(1), fly(1), 'r.') %lefts red
                       else
                           plot(flx(1), fly(1), '.', 'color', [1 .5 0]) %rights yellow
                       end
                   %end
                
               end
                %}
            
        end
    end
end

%erase unused pages (clusters) at the end of rate_matrix_4shuf
    for i = size(rate_matrix_4shuf,3):-1:1
        if nansum(nansum(rate_matrix_4shuf(:,:,i)))==0
            rate_matrix_4shuf(:,:,i) = [];
        end
    end


end
 