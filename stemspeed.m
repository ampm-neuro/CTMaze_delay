function [stem_speeds, indx] = stemspeed(eptrials, excluded)
%calculates how long it took the rat to run down the stem on every trial

%trials to output
trials = setdiff(unique(eptrials(:,6))', excluded);

stem_speeds = nan(length(trials), 1);


%iterate through trials
for i = 1:length(trials)
    
    trial = trials(i);
    
    stem_speeds(i,1:4) = [max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==1, 1)) - min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==1, 1))...
        max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==2, 1)) - min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==2, 1))...
        max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==3, 1)) - min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==3, 1))...
        max(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4, 1)) - min(eptrials(eptrials(:,6)==trial & eptrials(:,11)==4, 1))]; %section times
    stem_speeds(i,5) = trial; %trial number
    stem_speeds(i,6) = mode(eptrials(eptrials(:,6)==trial,9)); %accuracy
    
end

    %mean and std in each section
    mean_stm_sec = mean(stem_speeds(:,1:4));
    std_stm_sec = std(stem_speeds(:,1:4));
    
    %acceptably "ballistic"
    upper_lim = mean_stm_sec + std_stm_sec*3;
    
for i = 1:size(stem_speeds,1)
    
    if all((stem_speeds(i,1:4)-upper_lim)<0)
        stem_speeds(i, 7) = 1;
    else 
        stem_speeds(i, 7) = 0;
    end
    
end

%both correct AND ballistic
indx = sum(stem_speeds(:, [6 7]), 2)==2;
trials = trials(indx);

stem_speeds = nan(length(trials), 1);

%take average velocity
for i = 1:length(trials)
    
    trial = trials(i);
    
    stem_speeds(i,1) = mean(eptrials(~isnan(eptrials(:,12)) & eptrials(:,6)==trial & ismember(eptrials(:,11),1:4), 12))/100; %section velocity m/s
    
end




end