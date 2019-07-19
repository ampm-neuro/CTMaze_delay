function [safter, IDs] = ALL_keytimes_LR(stage_numbers, LR, evt)
%all performs some function copy pasted below on every session. typically
%for plotting


% stage numbers
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day

%set counters
%count = 0;

%min_pop = 8;
%popday = 0;

%prep
safter_start_hold = [];
safter_lighton_hold = [];
safter_stexit_hold = [];
safter_rwd_hold = [];
safter_start = [];
safter_lighton = [];
safter_stexit = [];
safter_rwd = [];

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);


%sesh_num = 0;

%iterate through subjects...and stages...and sessions
for subject = 1:length_subjects
    
    
    %figure;hold on
    
    %print update
    rat = file_names_subjects{:}(subject,1).name;
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
        length_sessions = size(file_list_sessions,1) -1;
    
        %iterate through sessions
        %for session = ceil(length_sessions/2):length_sessions
        for session = 1:length_sessions
            %session
            
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
            
            %
            if isempty(clusters) %|| length(clusters(:,1)) < min_pop
                continue
            else
                %popday = popday+1;
            end    

               %rejected
               all_trials = unique(eptrials(:,6));
               rejected_trials = reject_trials('all', lost_trials, out_stem, light_on_sect, [1 2.1 2.2 2.3 2.4], out_reward);
               included_trials = setdiff(all_trials,rejected_trials);
               error_trials = unique(eptrials(eptrials(:,9)==0,6));
               included_trials = setdiff(included_trials,error_trials);
               
               included_trials_L = setdiff(included_trials, unique(eptrials(eptrials(:,8)==2,6)));
               included_trials_R = setdiff(included_trials, unique(eptrials(eptrials(:,8)==1,6)));
               
               
               %remove one trial type
               if LR == 1
                   exclude_trl_type = unique(eptrials(eptrials(:,8)==2,6));
               elseif LR == 2
                   exclude_trl_type = unique(eptrials(eptrials(:,8)==1,6));
               end
               included_trials = setdiff(included_trials, exclude_trl_type);


               min_trials = 0;
               if length(included_trials_L) < min_trials || length(included_trials_R) < min_trials
                   continue
               end
               

               %for trial = included_trials((1+ceil(length(included_trials)/2)): end)'
               for trial = included_trials'
                   %event times
                   trial_start = min(eptrials(eptrials(:,6)==trial,1));
                   fl = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,13),[1 2]),1));
                   ent = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11),1),1));
                   ext = max(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,11),4),1));
                   rwd = min(eptrials(eptrials(:,6)==trial & ismember(eptrials(:,5),[1 2 3 4]),1));
                   
                   %key rates
                   if ~isempty(fl) && ~isempty(ent) && ~isempty(ext) && ~isempty(rwd)
                        safter_start_hold = [safter_start_hold ratewdw(eptrials, trial_start, trial_start+1, clusters)];
                        safter_lighton_hold = [safter_lighton_hold ratewdw(eptrials, fl, fl+1, clusters)];
                        safter_stexit_hold = [safter_stexit_hold ratewdw(eptrials, ext, ext+1, clusters)];
                        safter_rwd_hold = [safter_rwd_hold ratewdw(eptrials, rwd, rwd+1, clusters)];
                   end
               end
                
               %load
               %loads all cell events from the first trial up until the
               %trial in each session corresponding the the number of
               %elligible trials in the session with the least number of
               %elligible trials
                    
                    %conform trial number to minimum
                    if ~isempty(safter_start)
                        
                        if size(safter_start_hold,2) < size(safter_start,2)
                            safter_start = safter_start(:,1:size(safter_start_hold,2));
                            safter_lighton = safter_lighton(:,1:size(safter_start_hold,2));
                            safter_stexit = safter_stexit(:,1:size(safter_start_hold,2));
                            safter_rwd = safter_rwd(:,1:size(safter_start_hold,2));

                        elseif size(safter_start_hold,2) > size(safter_start,2)
                            safter_start_hold = safter_start_hold(:,1:size(safter_start,2));
                            safter_lighton_hold = safter_lighton_hold(:,1:size(safter_start,2));
                            safter_stexit_hold = safter_stexit_hold(:,1:size(safter_start,2));
                            safter_rwd_hold = safter_rwd_hold(:,1:size(safter_start,2));
                        end

                    end
                        
                   %load
                   safter_start = [safter_start; safter_start_hold];
                   safter_lighton = [safter_lighton; safter_lighton_hold];
                   safter_stexit = [safter_stexit; safter_stexit_hold];
                   safter_rwd = [safter_rwd; safter_rwd_hold];


                   %clear
                   safter_start_hold = [];
                   safter_lighton_hold = [];
                   safter_stexit_hold = [];
                   safter_rwd_hold = [];
                        

        end
    end
end

%{
safter_start_train = safter_start_train';
safter_lighton_train = safter_lighton_train';
%safter_stemtrance_train = safter_stemtrance_train';
safter_stexit_train = safter_stexit_train';
safter_rwd_train = safter_rwd_train';
safter_start_test = safter_start_test';
safter_lighton_test = safter_lighton_test';
%safter_stemtrance_test = safter_stemtrance_test';
safter_stexit_test = safter_stexit_test';
safter_rwd_test = safter_rwd_test';
%}
safter_start = safter_start';
safter_lighton = safter_lighton';
%safter_stemtrance = [safter_stemtrance ratewdw(eptrials, ent, ent+1, clusters)];
safter_stexit = safter_stexit';
safter_rwd = safter_rwd';



%condense outputs
if evt == 0
    safter = [safter_start; safter_lighton; safter_stexit; safter_rwd]; %safter_stemtrance_train; 
    IDs = [ones(size(safter_start,1), 1); repmat(2, size(safter_lighton, 1), 1); repmat(3, size(safter_stexit, 1), 1); repmat(4, size(safter_rwd, 1), 1)];

elseif evt == 1
    safter = safter_start; %safter_stemtrance_train; 
    IDs = ones(size(safter_start,1), 1);

elseif evt == 2
   	safter = safter_lighton; %safter_stemtrance_train; 
    IDs = ones(size(safter_lighton,1), 1);

elseif evt ==3
  	safter = safter_stexit; %safter_stemtrance_train; 
    IDs = ones(size(safter_stexit,1), 1);

elseif evt == 4
 	safter = safter_rwd; %safter_stemtrance_train; 
    IDs = ones(size(safter_rwd,1), 1);

end





%train_IDs = [ones(size(safter_lighton_train,1), 1)];
%test_IDs = train_IDs;

%classify subsets and average (to deal with covariance estimation limits)
%class1 = classify(safter_test(:, 1:50), safter_train(:, 1:75), train_IDs);
%class2 = classify(safter_test(:, 51:100), safter_train(:, 76:150), train_IDs);
%class3 = classify(safter_test(:, 101:150), safter_train(:, 151:225), train_IDs);
%class4 = classify(safter_test(:, 151:200), safter_train(:, 151:200), train_IDs);

%class = classify(safter_test, safter_train, train_IDs, 'diaglinear');

%means
%prop_correct = [sum(class1==test_IDs) sum(class2==test_IDs) sum(class3==test_IDs) sum(class4==test_IDs) sum(class5==test_IDs)]./[length(class1) length(class2) length(class3) length(class4) length(class5)];
%prop_correct = class;

end


function [rate] = ratewdw(eptrials, low, high, clusts)
%calculates rates for each cell in clusts during the time window from low
%to high

    rate = histc(eptrials(eptrials(:,1)>=low & eptrials(:,1)<=high, 4), clusts);
    rate = rate./repmat(high-low, size(rate));

end
 