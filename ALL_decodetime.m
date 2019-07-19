function [rate_matrices, posterior_means, posterior_trials, reshape_trials, bin_times, rm_bins_normal_means, velocity_vect, acceleration_vect] = ALL_decodetime(stage_numberss, wndw, bins, boxcar_shift) %, posterior_means, posterior_trials, reshape_trials
%all performs some function copy pasted below on every session. typically
%for plotting


% stage numbers
%   1 = crit/overtrain days
%   2 = middle learning day
%   3 = first learning day


%start [.1 5.5700]
%light [-1.7312 3.5599]
%entrance [-2.0734 3.1905]
%exit [-2.8278 2.4359]
%reward [-4.2911 1]

%set counters
%count = 0;
rate_matrices = cell(26,1);
velocity_vect = [];
acceleration_vect = [];
sesh_count = 0;

bin_times = boxcar_bintimes(wndw, bins, boxcar_shift);

%cramped coding of folder access. see ALL.m file types for documentation
file_list_subjects=dir('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/');
file_list_subjects(1:2) = [];
file_names_subjects={file_list_subjects([file_list_subjects(:).isdir])};
length_subjects = size(file_names_subjects{:},1);

sesh_num = 0;

%iterate through subjects...and stages...and sessions
for subject = 1:length_subjects
    %print update
    rat = file_names_subjects{:}(subject,1).name;
    %get all the things in subject folder...
    file_list_stages = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name)));
    file_list_stages(1:2) = [];
    file_names_stages = {file_list_stages([file_list_stages(:).isdir])};
    length_stages = size(file_names_stages{:},1);
    
    %stage_numbers
    
    %iterate through stages...and sessions
    for stage = stage_numberss
        
        %print update
        task = file_names_stages{:}(stage,1).name;
        %get all the *.mat in subject folder...
        file_list_sessions = dir(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/', num2str(file_names_subjects{:}(subject,1).name), '/', num2str(file_names_stages{:}(stage,1).name, '/*.mat')));
        file_list_sessions(1:2) = [];
        length_sessions = size(file_list_sessions,1);
    
        %iterate through sessions
        for session = 1:length_sessions
            rate_mtx = [];
            
            if session <10
                session_num = strcat('0', num2str(session));
            else
                session_num = session;
            end
            
            %load
            eptrials = [];
            origin_file = [];
            clusters = [];
            light_on_sect = [];
            load(strcat('/Users/ampm/Documents/MATLAB/lindseyvedder/neurodata/',num2str(rat),'/' ,num2str(task),'/' ,num2str(session_num), '.mat'));
            
            if ~isempty(clusters)
                
               
                %correct trials
                corrects = unique(eptrials(eptrials(:,9)==1,6));
                
                %flag times
                %flag_events = nan(size(corrects));
                for trl_num = 1:length(corrects)

                    %event time
                    %et = min(eptrials(eptrials(:,6)==trl_num, 1)); %trial start
                    et = min(eptrials(eptrials(:,6)==trl_num & ismember(eptrials(:,13), [1 2]), 1)); %lighton
                    %et = min(eptrials(eptrials(:,6)==trl_num & ismember(eptrials(:,11), 1), 1)); %stem entrance
                    %et = max(eptrials(eptrials(:,6)==trl_num & ismember(eptrials(:,11), 4), 1)); %stem exit
                    %et = min(eptrials(eptrials(:,6)==trl_num & ismember(eptrials(:,5), [1 2 3 4]), 1)); %reward

                    if ~isempty(et)
                        %for now, light on
                        %flag_events(trl_num) = et;
                        [rate_mtx_part] = spike_count_boxcar(eptrials(:, [1 4]), clusters, bin_times, et);
                        rate_mtx = cat(3, rate_mtx, rate_mtx_part);
                        
                        [accel_mtx_part] = spike_count_boxcar(eptrials(:, [1 14]), 1, bin_times, et);
                        [veloc_mtx_part] = spike_count_boxcar(eptrials(:, [1 12]), 1, bin_times, et);
                        velocity_vect = cat(3, velocity_vect, veloc_mtx_part);
                        acceleration_vect = cat(3, acceleration_vect, accel_mtx_part);
                        
                    end
                end

                sesh_count = sesh_count + 1;
                rate_matrices(sesh_count) = {rate_mtx};
                
                
            end
        end
    end
end



%merge rate_matrices
[rm_bin_means_train, rm_bins_test] = merge_rate_matrices(rate_matrices);

    %create normalized copy for later plot
    rm_bins_normal = reshape(rm_bins_test, size(rm_bins_test,1), size(rm_bins_test,2)*size(rm_bins_test,3));
    cell_means = mean(rm_bins_normal,2);
    cell_std = std(rm_bins_normal, [], 2);

    rm_bins_normal = rm_bins_normal - repmat(cell_means, 1, size(rm_bins_normal, 2));
    rm_bins_normal = rm_bins_normal ./ repmat(cell_std, 1, size(rm_bins_normal, 2));

    rm_bins_normal = reshape(rm_bins_normal, size(rm_bins_test,1), size(rm_bins_test,2), size(rm_bins_test,3));
    rm_bins_normal_means = mean(rm_bins_normal,3);
    
    mean_normal_rates = mean(rm_bins_normal_means);
    std_normal_rates = std(mean_normal_rates);

IDs = repmat(1:size(rm_bin_means_train,2), 1, 15);
training = reshape(rm_bin_means_train, [246 , size(rm_bins_test,2)*15]);
[posterior_trials] = decode_time(training, rm_bins_test, IDs, bins);

count = 0; 
count2 = 1; 
for i = 1:20
    reshape_trials(1:size(rm_bins_test,2), :, count2) = posterior_trials(count+1:count+size(rm_bins_test,2), :);
    count = count+size(rm_bins_test,2); 
    count2 = count2+1; 
end

%means and stds
posterior_means = mean(reshape_trials,3);
counts_means = mean(rm_bins_normal,3);
posterior_means_diag = []; 
counts_means_diag = [];
for i = 1:size(reshape_trials,1)
    posterior_means_diag = [posterior_means_diag; std(squeeze(reshape_trials(i,i,:)))];
    counts_means_diag = [counts_means_diag; std(squeeze(rm_bins_normal(i,i,:)))];
end

    
%PLOTS
%
%median event BIN times, based on actual times
%
%trial start
%{
light_on = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 1.8312);
stem_ent = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 2.2383);
stem_ext = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 3.0029);
rwd_time = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 4.5700);
%}

%light on
%
light_on = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0);
stem_ent = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0.3422);
stem_ext = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 1.0966);
rwd_time = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 2.5599);
%}

%stem entrance
%{
light_on = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -0.3422);
stem_ent = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0);
stem_ext = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0.7514);
rwd_time = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 2.1905);
%}

%stem exit
%{
light_on = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -1.0966);
stem_ent = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -0.7514);
stem_ext = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0);
rwd_time = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 1.4359);
%}

%reward
%{
light_on = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -2.5599);
stem_ent = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -2.1905);
stem_ext = interp1(bin_times(:,2), 1:length(bin_times(:,2)), -1.4359);
rwd_time = interp1(bin_times(:,2), 1:length(bin_times(:,2)), 0);
%}



%axes
tickpoints = 1:((size(posterior_means,1)-1)/10):size(posterior_means,1);
ticklabels = linspace(bin_times(1,2), bin_times(end,2), length(tickpoints));

%probability matrix (heatmap)
figure; hold on
imagesc(posterior_means'); colorbar; colormap jet; 
caxis([0 .4])
plot([light_on light_on], [0.5 size(posterior_means,1)+.5], 'r-')
plot([stem_ent stem_ent], [0.5 size(posterior_means,1)+.5], 'k-')
plot([stem_ext stem_ext], [0.5 size(posterior_means,1)+.5], 'k-')
plot([rwd_time rwd_time], [0.5 size(posterior_means,1)+.5], 'k-')

ax = gca;
ax.XTick = tickpoints;
ax.XTickLabel = round(ticklabels,3);
ax.YTick = tickpoints;
ax.YTickLabel = round(ticklabels,3);
axis square
axis([.5 size(posterior_means,1)+.5 .5 size(posterior_means,1)+.5])
set(gca,'TickLength',[0, 0]);

%diaganol plot
figure; hold on

plot(smooth(posterior_means(diag_mask(size(posterior_means,1))),3)); title('smoothed')
plot(smooth(posterior_means(diag_mask(size(posterior_means,1))) - posterior_means_diag./sqrt(20),3), 'b');
plot(smooth(posterior_means(diag_mask(size(posterior_means,1))) + posterior_means_diag./sqrt(20),3), 'b');

xlim([0 size(posterior_means,1)+1])

plot([light_on light_on], [0 1], 'r-')
plot([stem_ent stem_ent], [0 1], 'k-')
plot([stem_ext stem_ext], [0 1], 'k-')
plot([rwd_time rwd_time], [0 1], 'k-')

ax = gca;
ax.XTick = tickpoints;
ax.XTickLabel = round(ticklabels,3);



%normalized rates plot
figure; hold on

plot(smooth(mean_normal_rates,3)); title('smoothed')
plot(smooth(mean_normal_rates - std_normal_rates./sqrt(size(std_normal_rates,1)),3), 'b');
plot(smooth(mean_normal_rates + std_normal_rates./sqrt(size(std_normal_rates,1)),3), 'b');

xlim([0 size(posterior_means,1)+1])
plot([0 size(posterior_means,1)+1], [0 0], 'k--')

plot([light_on light_on], [-1 1], 'r-')
plot([stem_ent stem_ent], [-1 1], 'k-')
plot([stem_ext stem_ext], [-1 1], 'k-')
plot([rwd_time rwd_time], [-1 1], 'k-')

ax = gca;
ax.XTick = tickpoints;
ax.XTickLabel = round(ticklabels,3);


%normalized INDIVIDUAL rates plot
figure; hold on

subsamp = randperm(size(rm_bins_normal,1));
plot(rm_bins_normal_means(subsamp(1:50),:)'); title('smoothed')
%plot(smooth(mean_normal_rates - std_normal_rates./sqrt(size(std_normal_rates,1)),3), 'b');
%plot(smooth(mean_normal_rates + std_normal_rates./sqrt(size(std_normal_rates,1)),3), 'b');

xlim([0 size(posterior_means,1)+1])
plot([0 size(posterior_means,1)+1], [0 0], 'k--')

plot([light_on light_on], [-3 3], 'r-')
plot([stem_ent stem_ent], [-3 3], 'k-')
plot([stem_ext stem_ext], [-3 3], 'k-')
plot([rwd_time rwd_time], [-3 3], 'k-')

ax = gca;
ax.XTick = tickpoints;
ax.XTickLabel = round(ticklabels,3);

%acceleration / velocity plot
%
figure; hold on

plot(mean(velocity_vect,3), 'r'); title('accel and veloc')
plot(mean(velocity_vect,3) - std(velocity_vect, [], 3)./sqrt(size(velocity_vect,3)), 'r');
plot(mean(velocity_vect,3) + std(velocity_vect, [], 3)./sqrt(size(velocity_vect,3)), 'r');

plot(mean(acceleration_vect,3), 'y')
plot(mean(acceleration_vect,3) - std(acceleration_vect, [], 3)./sqrt(size(acceleration_vect,3)), 'y');
plot(mean(acceleration_vect,3) + std(acceleration_vect, [], 3)./sqrt(size(acceleration_vect,3)), 'y');

xlim([0 size(posterior_means,1)+1])
plot([0 size(posterior_means,1)+1], [0 0], 'k--')

plot([light_on light_on], [-75 75], 'r-')
plot([stem_ent stem_ent], [-75 75], 'k-')
plot([stem_ext stem_ext], [-75 75], 'k-')
plot([rwd_time rwd_time], [-75 75], 'k-')

ax = gca;
ax.XTick = tickpoints;
ax.XTickLabel = round(ticklabels,3);
%}



end
 