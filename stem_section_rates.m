function ssrs = stem_section_rates(eptrials, clusters, rejected_trials, idx)
%calculates the firing rate of each cell in clusters in each section of the
%stem on each trial. Outputs a 3d matrix with rows for trials,
%columns for stem section, and pages for cell. Output can be used to
%calculate Emma Wood splitters.
%
%ignores trials in rejected_trials vector
%


%accepted trials
%trials = setdiff(1:max(eptrials(:,6)), rejected_trials);
trials = 1:max(eptrials(:,6));
trials = trials(idx);

%number of sections
num_stem_sec = max(eptrials(~isnan(eptrials(:,11)),11));

%preallocate output (trials, stem sections, clusters)
ssrs = nan(length(trials), num_stem_sec, length(clusters));

for row = 1:length(trials)
    trial = trials(row);
    
    %find first time point of each section, and last fourth (exit)
    first = min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==1,1));
    second = min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==2,1));
    third = min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==3,1));
    fourth = min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4,1));
    exit = max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4,1));
    
    if isempty(first)
        continue
    end
    
    %time indices
    t1_idx = eptrials(:,1)>=first & eptrials(:,1)<second;
    t2_idx = eptrials(:,1)>=second & eptrials(:,1)<third;
    t3_idx = eptrials(:,1)>=third & eptrials(:,1)<fourth;
    t4_idx = eptrials(:,1)>=fourth & eptrials(:,1)<=exit;
    
    %find time (in seconds) spent in each section
    t1 = length(eptrials(eptrials(:,4)==1 & t1_idx, 4))*.01;
    t2 = length(eptrials(eptrials(:,4)==1 & t2_idx, 4))*.01;
    t3 = length(eptrials(eptrials(:,4)==1 & t3_idx, 4))*.01;
    t4 = length(eptrials(eptrials(:,4)==1 & t4_idx, 4))*.01;
    
    %count spikes in each section 
    for page = 1:length(clusters)
        cell = clusters(page);

        s1 = length(eptrials(eptrials(:,4)==cell & t1_idx, 4));
        s2 = length(eptrials(eptrials(:,4)==cell & t2_idx, 4));
        s3 = length(eptrials(eptrials(:,4)==cell & t3_idx, 4));
        s4 = length(eptrials(eptrials(:,4)==cell & t4_idx, 4));
        
        %load rates into output matrix
        ssrs(row, :, page) = [s1/t1 s2/t2 s3/t3 s4/t4];
        
    end
end