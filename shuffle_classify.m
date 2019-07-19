function prop_correct = shuffle_classify(num_shufs)
%repeatedly runs ALL_keytimes

prop_correct = [];

for i = 1:num_shufs
    
    i
    
    [~, ~, ~, ~, prop_correct_part] = ALL_keytimes(1);
    prop_correct = [prop_correct; prop_correct_part];
    
end



end