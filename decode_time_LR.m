function [posterior_trials] = decode_time_LR(rm_bins_training, rm_bins_test, IDs, bins)
%calculate the probability that the rat is currently occupying the 0-bin_lengths bin
%after light onset, from among the BINS bins before and BINS bins after light
%onset.
    posterior_trials = [];
    
    %for each trial
    for trial = 1:size(rm_bins_test,3)

        test = rm_bins_test(:,:,trial);
        test = test';
       
        
        %classify
        [~, ~, posterior] = classify(test, rm_bins_training', IDs', 'diagLinear');
        %
        %[~, posterior] = bayesian_decode(test, rm_bins_training', IDs', bins);
        
        posterior_trials = [posterior_trials; posterior];
        
    end

end