function sd = standard_distance(trials1, trials2)
%standard_distance calculates the difference between the mean of trials1
%and the mean of trials2 divided by the pooled standard deviation

m1 = mean(trials1);
m2 = mean(trials2);

dif1 = trials1 - repmat(m1,size(trials1));
dif2 = trials2 - repmat(m2,size(trials2));

pstd = mean([dif1; dif2])/sqrt(length([trials1;trials2]));


sd = abs(m1-m2)/pstd;

end