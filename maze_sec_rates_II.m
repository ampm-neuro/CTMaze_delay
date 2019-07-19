function folded_sec_rates = maze_sec_rates_II(eptrials, cluster, varargin)
%another attempt to measure the rates observed in each section on each
%trial


%a list of rejected trials isnt required, but allowed
if ~isempty(varargin)
    rejected_trials =  varargin{1};
else
    rejected_trials = [];
end


%preallocate ([folded rates, choice, accuracy])
folded_sec_rates = nan(max(eptrials(:,6)), 7);

%accepted trials
trials = setdiff(1:max(eptrials(:,6)), rejected_trials);

%iterate through trials
for trial = trials
    
    %calculate spikes and times
    sec_spikes = histc(eptrials(eptrials(:,6)==trial & eptrials(:,4)==cluster, 10), 1:7);
        if isrow(sec_spikes)
            sec_spikes = sec_spikes';
        end    
    sec_times = histc(eptrials(eptrials(:,6)==trial & eptrials(:,4)==1, 10), 1:7).*.01;
    
    %calculate rates (spikes/times)
    unfolded_sec_rates = sec_spikes./sec_times;
    
    %load rates and trial choice into folded sections
    if mode(eptrials(eptrials(:,6)==trial, 8)) == 1
        folded_sec_rates(trial, 1:5) = [unfolded_sec_rates(1:3)' unfolded_sec_rates(4) unfolded_sec_rates(6)];
        folded_sec_rates(trial, 6) = 1;
    elseif mode(eptrials(eptrials(:,6)==trial, 8)) == 2
        folded_sec_rates(trial, 1:5) = [unfolded_sec_rates(1:3)' unfolded_sec_rates(5) unfolded_sec_rates(7)];
        folded_sec_rates(trial, 6) = 2;
    end
    
    %load accuracy
    if mode(eptrials(eptrials(:,6)==trial, 9)) == 0
        folded_sec_rates(trial, 7) = 0;
    elseif mode(eptrials(eptrials(:,6)==trial, 9)) == 1
        folded_sec_rates(trial, 7) = 1;
    end
end
