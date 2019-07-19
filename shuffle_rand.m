function out = shuffle_rand(trials, shuffs)

test = [zeros(trials/2,1); ones(trials/2, 1)];
out = nan(shuffs, 1);

for i = 1:shuffs

    shuf = test(randperm(length(test)));
    
    out(i) = sum(shuf==test);
end

out = sort(out);


idx_prep = 1:length(out);



twoptfive_pct = find(abs(idx_prep - repmat(.025*length(idx_prep), 1, length(idx_prep))) == min(abs(idx_prep - repmat(.025*length(idx_prep), 1, length(idx_prep)))), 1, 'first')
low = out(twoptfive_pct)/trials

mean_out = mean(out)/trials

ninetysevenptfive_pct = find(abs(idx_prep - repmat(.975*length(idx_prep), 1, length(idx_prep))) == min(abs(idx_prep - repmat(.975*length(idx_prep), 1, length(idx_prep)))), 1, 'first')
high = out(ninetysevenptfive_pct)/trials

end