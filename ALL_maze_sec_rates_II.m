function [rate_output LR_output] = ALL_maze_sec_rates_II(stage_numbers)
%iterates through every session file and cluster, calculating trial-by-trial
%maze section firing rates
%
%builds a (accurate trials X 5 folded sections X number of clusters)
%3d matrix of folded firing rates for each session with a page for
%each cluster and only rows for correct trials
%
%adds each of these session matrices to one overall cell matrix
%
%a second cell matrix holds the L's and R's of correct trials corresponding
%to each cell in the overall cell matrix

% stage number input options
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day

%set counters
sesh_num = 0;
cell_count = 0;

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);


%iterate through subjects...and stages...and sessions...and clusters
for subject = 1:length_subjects
    %print update
    rat = file_names_subjects{:}(subject,1).name
    %get all the things in subject folder...
    file_list_stages = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name)));
    file_list_stages(1:2) = [];
    file_names_stages = {file_list_stages([file_list_stages(:).isdir])};
    length_stages = size(file_names_stages{:},1);
    
    %iterate through stages...and sessions...and clusters
    for stage = stage_numbers
        %print update
        task = file_names_stages{:}(stage,1).name
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(file_names_stages{:}(stage,1).name, '/*.mat')));
        file_list_sessions(1:2) = [];
        length_sessions = size(file_list_sessions,1);
    
        %iterate through sessions...and clusters
        for session = 1:length_sessions
            
            %for loading session file
            if session <10
                session_num = strcat('0', num2str(session))
            else
                session_num = session
            end
            
            %initialize session variables
            eptrials = [];
            origin_file = [];
            clusters = [];
            light_on_sect = [];
            lost_trials = [];
            
            %load session file
            load(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/',num2str(rat),'/' ,num2str(task),'/' ,num2str(session_num), '.mat'));
            
            %skip cluster-less sessions
            if ~isempty(clusters)
                
                %keep track of number of unique sessions
                sesh_num = sesh_num+1;

                %preallocate rate_matrix for session
                %num_correct_trials = length(unique(eptrials(eptrials(:,9)==1, 6)));
                num_correct_trials = length(unique(eptrials(:, 6)));
                session_rates = nan(num_correct_trials - length(lost_trials), 5, length(clusters));
                
                %iterate through clusters
                for cell_num = 1:length(clusters)
                    clust = clusters(cell_num);
                    
                    %keep track of number of unique cells
                    cell_count = cell_count + 1;
                                    
                    %calculate folded rates, choices, accuracy
                    folded_sec_rates = maze_sec_rates_II(eptrials, clust, lost_trials);
                    
                    %pull out accurate-only LRs for shuffle
                    if cell_num == 1
                        %session_LRs{sesh_num} = folded_sec_rates(folded_sec_rates(:,7)==1,6);
                        session_LRs{sesh_num} = folded_sec_rates(ismember(folded_sec_rates(:,7), [0 1]),6);
                    end
                    
                    %load matrix of correct-trial rates for all cells
                    %session_rates(:,:,cell_num) = folded_sec_rates(folded_sec_rates(:,7)==1,1:5);
                    session_rates(:,:,cell_num) = folded_sec_rates(ismember(folded_sec_rates(:,7), [0 1]),1:5);
                end
                
                
                %add cell with session (cluster) rates to rate matrix
                rate_matrix_II{sesh_num} = session_rates;
                
            end
        end
        
        rate_output = rate_matrix_II;
        LR_output = session_LRs;
        
        %overall_output{stage} = stage_output;
        
    end

end

end