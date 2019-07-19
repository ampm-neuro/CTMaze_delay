function [msrs msrs_fold] = maze_section_rates(eptrials, clusters, lost_trials, light_on_sect, varargin)
%calculates the standard distance between firing rates in corresponding
%maze sections of the left and right trial types. ONLY INCLUDES TIMES WHEN
%THE LIGHT IS ON

msrs_fold = [];

%variable arguments in include the accepted light-ons and the session
%number
if ~isempty(varargin)
   accept_lights = varargin{1}; 
    %accepted trials - nix lost trials and trials where the light doesnt turn on
    rejected_trials = unique([lost_trials; light_on_sect(not(ismember(light_on_sect(:,2),accept_lights)),1)]);
    %
    %trial choice types
    class_mtx = trl_class(eptrials, rejected_trials);
 
    if nargin == 6      
        session_number = varargin{2};
    else
        session_number = 1;
    end
end

%accepted trials - nix lost trials and trials where the light doesnt turn on
if exist('rejected_trials', 'var')
    trials = setdiff(1:max(eptrials(:,6)), rejected_trials);
else
    trials = setdiff(1:max(eptrials(:,6)), unique([lost_trials; light_on_sect(light_on_sect(:,2)==0,1)]));
end

%number of sections
num_maze_sec = 7;
num_folded_sec = 5;

%preallocate output (trials, stem sections, clusters)
msrs = nan(max(eptrials(:,6)), num_maze_sec, length(clusters));
msrs_fold = nan(max(eptrials(:,6)), num_folded_sec+4, size(msrs,3));

%iterate through trials
for row = trials
    
    %set desired time range
        %first instant of light-on
        first = min(eptrials(eptrials(:,6)==row & ~isnan(eptrials(:,13)),1));
        %trial_start = min(eptrials(eptrials(:,6)==row & eptrials(:,4)==1,1));
        %first = trial_start;

        %light off or rwd+3s, whichever comes first
        last_time = max(eptrials(eptrials(:,6)==row & ismember(eptrials(:,5), 1:4) ,1)) + 3;
        %last_light = max(eptrials(eptrials(:,6)==row & ~isnan(eptrials(:,13)),1));
        %last = min([last_time last_light]);
        last = last_time;

        %time index
        t_idx = eptrials(:,1)>=first & eptrials(:,1)<=last;
    
    %find time (in seconds) spent in each section on that trial (w/in time
    %range). will be used for all clusters
    time = nan(num_maze_sec,1);
    for sect = 1:num_maze_sec
        time(sect) = length(eptrials(eptrials(:,6)==row & eptrials(:,10)==sect & eptrials(:,4)==1 & t_idx, 4))*.01;
    end
    
    %nans for non-visted sections
    time(time==0) = nan;
        
    %count spikes and divide by time in each section on that trial
    for page = 1:length(clusters)
        cell = clusters(page);
        
        for sect = 1:num_maze_sec
            %load rates into output matrix
            msrs(row, sect, page) = length(eptrials(eptrials(:,6)==row & eptrials(:,10)==sect & eptrials(:,4)==cell & t_idx, 4))/time(sect);
            %msrs(trial, section, cluster)
        end        
        
        if ~isempty(varargin)
            %5 column rate matrix (folded)
            msrs_fold(1:size(msrs,1),1:3, page) = msrs(:, 1:3, page);
            msrs_fold(~isnan(msrs(:,4)),4, page) = msrs(~isnan(msrs(:,4)),4, page);
            msrs_fold(~isnan(msrs(:,5)),4, page) = msrs(~isnan(msrs(:,5)),5, page);
            msrs_fold(~isnan(msrs(:,6)),5, page) = msrs(~isnan(msrs(:,6)),6, page);
            msrs_fold(~isnan(msrs(:,7)),5, page) = msrs(~isnan(msrs(:,7)),7, page);
            msrs_fold(1:size(msrs,1),6, page) = class_mtx(:,3); %choice
            msrs_fold(1:size(msrs,1),7, page) = page; %cell number (not id)
            msrs_fold(1:size(msrs,1),8, page) = session_number; %session
            msrs_fold(1:size(msrs,1),9, page) = class_mtx(:,4); %accuracy
        end
        
    end
end

end