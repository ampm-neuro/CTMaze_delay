function [rm_bin_means_train, rm_bins_test] = merge_rate_matrices(rate_matrices)

rm_bin_means_train = cell(size(rate_matrices));

    for sesh = 1:length(rate_matrices)

        %train_sesh_rng = 1:(size(rate_matrices{sesh}, 3)-20);
        %rm_bin_means_train(sesh) = {mean(rate_matrices{sesh}(:,:,train_sesh_rng),3)};
        
        train_sesh_rng = 1:15;
        rm_bin_means_train(sesh) = {rate_matrices{sesh}(:,:,train_sesh_rng)};

    end
    
rm_bin_means_train = cell2mat(rm_bin_means_train);
    
    
rm_bins_test = cell(size(rate_matrices));

    for sesh = 1:length(rate_matrices)

        test_sesh_rng = (size(rate_matrices{sesh}, 3)-19):size(rate_matrices{sesh}, 3);
        
        rm_bins_test(sesh) = {rate_matrices{sesh}(:,:,test_sesh_rng)};

    end

rm_bins_test = cell2mat(rm_bins_test); 
    

end