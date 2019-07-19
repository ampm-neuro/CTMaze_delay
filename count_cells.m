function [cell_counts, cell_type_counts, type_combos, standard_distances, light_counts, interaction_counts, sect_dist_rank, rwd_count, light_count, split_pvals, CA_pvals] = count_cells(min_pop_size, stage_numbers, accept_light)
% count_cells calculates and counts the number of cells of different types
%
% It requires a specific data organization structure
%
% It only consideres clusters with at least confidence_level confidence,
% and populations of at least size min_pop_size. It also only looks at
% learning stages indicated by stage_numbers, where: 
%
% stage numbers
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day
%   4 = pretraining day
%

%set counters
cell_counts = [0;0]; %cells; populations
cell_type_counts = [0;0;0;0;0]; %rwd; place; stem_diff; choice_diff; arm_diff
type_combos = [0;0;0]; %rwd-place; rwd-stem; place-stem
standard_distances = nan(360, 5);
light_counts = [];
interaction_counts = [0;0;0]; %interaction_only_count; interaction_and_main_count; main_only_count
sect_dist_rank = [];
rwd_count = [0;0];
light_count = [0;0];
split_pvals = [];
CA_pvals = [];


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
        task = file_names_stages{:}(stage,1).name
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(file_names_stages{:}(stage,1).name, '/*.mat')));
        file_list_sessions(1:2) = [];
        length_sessions = size(file_list_sessions,1);
    
        %iterate through sessions
        for session = 1:length_sessions
            
            if session <10
                session_num = strcat('0', num2str(session))
            else
                session_num = session
            end
            
            
            eptrials = [];
            origin_file = [];
            clusters = [];
            light_on_sect = [];
            load(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/',num2str(rat),'/' ,num2str(task),'/' ,num2str(session_num), '.mat'));
                        
            %screen out cell-less sessions
            if ~isempty(clusters)
                
                %concantonate light counts
                light_counts = [light_counts; light_on_sect(:,2)];
                
                
                %calculate standard distances
                %std_dists = sect_dists(eptrials, clusters, lost_trials, light_on_sect, accept_light);
                %standard_distances(cell_counts(1)+1:cell_counts(1)+length(clusters),:) = std_dists;
                
                %number of cells
                cell_counts(1) = cell_counts(1) + length(clusters);
                
                %number of rwd discriminating cells
                summary = light_rwd_cell(eptrials, clusters, .5, .5, 0);
                rwd_count(1) = rwd_count(1) + nansum(summary(:,3)); %rwd_encode_count
                rwd_count(2) = rwd_count(2) + nansum(summary(:,4)); %rwd_discrim_count
                
                %number of light discriminating cells
                light_count(1) = light_count(1) + nansum(summary(:,1)); %light_encode_count
                light_count(2) = light_count(2) + nansum(summary(:,2)); %light_discrim_count
                
                    
                %number of place cells
                %dms_place = spatialfield_batch(eptrials, clusters);
                %place_cell_count  = place_cell_count + sum(dms_place);
                
                %number of stem-differenting (last input is the acceptable
                %areas for the light to turn on. 1 is start (standard
                %protocol), 2.1 is the first stem area, 2.2 is the second
                %and so on.
                [stem_diff, interaction_only_count, interaction_and_main_count, main_only_count, ~, ~, split_pval_lookup] = splitter_test(eptrials, clusters, lost_trials, out_stem, light_on_sect, accept_light, subject, session);
                cell_type_counts(3) = cell_type_counts(3) + sum(stem_diff);
                interaction_counts(1:3) = interaction_counts(1:3) + [sum(interaction_only_count); sum(interaction_and_main_count); sum(main_only_count)];
                %sect_dist_rank = [sect_dist_rank; ranks];
                split_pvals = [split_pvals; split_pval_lookup];
                
                
                
                %number of chc/arm/rwd differentiating cells
                [CA_diff CA_pval_lookup] = arm_anova(eptrials, clusters, subject, session);
                cell_type_counts(4) = cell_type_counts(4) + nansum(CA_diff(:,1));
                cell_type_counts(5) = cell_type_counts(5) + nansum(CA_diff(:,2));
                CA_pvals = [CA_pvals; CA_pval_lookup];

                
                %reward cells that are also differentiating on stem
                %rwd_stem_comb_rwdcell = rwd_stem_comb_rwdcell + sum(sum([rwd_cell stem_diff],2)>1);
                
                %place cells that are also differentiating on stem
                %place_stem_comb = place_stem_comb + sum(sum([dms_place stem_diff],2)>1);
                
                %count populations
                if length(clusters)>=min_pop_size %CONTROL POPULATION SIZE
                    cell_counts(2) = cell_counts(2) + 1;
                end
            end
        end
    end
end
end






                
                
                
              