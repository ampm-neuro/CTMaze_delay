function class_mtx = trl_class(eptrials, idx)%varargin)
%outputs a matrix with trial numbers in the first column and a trial-
%by-trial classification in 3 more columns for: type (1/2), choice (1/2),
%accuracy (0/1)

    trials = 1:max(eptrials(:,6));
    trials = trials(idx);

    %preallocate
    class_mtx = nan(max(eptrials(:,6)),4);
    
    %set trial numbers
    class_mtx(:,1) = 1:max(eptrials(:,6));
    
    %if nargin == 2
        %rejected_trials = varargin{1};
    %    idx = varargin{1};
    %else
    %    rejected_trials = [];
    %end

    for trial = trials%setdiff(1:max(eptrials(:,6)), rejected_trials)
        class_mtx(trial, 2) = mode(eptrials(eptrials(:,6)==trial, 7));%type
        class_mtx(trial, 3) = mode(eptrials(eptrials(:,6)==trial, 8));%choice
        class_mtx(trial, 4) = mode(eptrials(eptrials(:,6)==trial, 9));%accuracy
    end
        
    %lr_mtx(isnan(lr_mtx(:,2)), :)=[];

end