
function velocity = plottrials_veloc(varargin)
%function plottrials(eptrials, cell_range)
%plot the trajectory for all trials, or just a specified range

figures = 1;
high= 35;
low = 25;


%check inputs
switch nargin
    case 0
        error(message('need more inputs'))
    case 1
        eptrials = varargin{1};
    case 2
        eptrials = varargin{1};
        trials = varargin{2};   
    otherwise
        error(message('too many inputs'))
end

%what trials to plot
if exist('trials', 'var')
    range = trials;
else
    range = 1:1:(max(eptrials(:,6)));
end


%establish distance units
rwd_dist = pdist([mean(eptrials(eptrials(:,5)==1,2)) mean(eptrials(eptrials(:,5)==1,3)); mean(eptrials(eptrials(:,5)==2,2)) mean(eptrials(eptrials(:,5)==2,3))]);
cm_per_matlab_unit = 115.5/rwd_dist;

%calulate velocity
for trial = range
    
    local_eptrials = eptrials(eptrials(:,6)==trial & eptrials(:,4)==1, :);

    veloc = zeros(length(local_eptrials(:,1)-1),1);
    
    p1 = local_eptrials(1:end-1, 2:3);
    p2 = local_eptrials(2:end, 2:3);
    t1 = local_eptrials(1:end-1, 1);
    t2 = local_eptrials(2:end, 1);
    
    
    for instant = 1:length(local_eptrials(:,1))-1

            veloc(instant+1) = pdist([p1(instant,:); p2(instant,:)])/(t2(instant)-t1(instant));
            
    end

    veloc = veloc.*cm_per_matlab_unit;
    time = local_eptrials(:,1) - repmat(local_eptrials(1,1), size(local_eptrials(:,1)));
    
    %smooth
    veloc = smooth(veloc, 50);
    
    %output
    velocity = [local_eptrials(:,1) veloc];
    
    
    %earliest rwd departure
    rwd_time_allotment = 3; %arbitrary
    post_rwd = eptrials(eptrials(:,5)>0 & eptrials(:,6)==trial,1) + rwd_time_allotment  - local_eptrials(1,1);
    %end of iti
    last_iti = max(eptrials(eptrials(:,10)==8 & eptrials(:,6)==trial,1)) - local_eptrials(1,1);

    
    if figures == 1
        figure
        plot(time(2:end), veloc(2:end), 'k-', 'LineWidth',1)
        hold on
        plot([time(2) time(end)], [high high],'r-', 'LineWidth',1)
        plot([time(2) time(end)], [low low],'r-', 'LineWidth',1)
        plot([post_rwd post_rwd], [0 max(veloc(2:end))*1.2],'b-', 'LineWidth',1)
        plot([last_iti last_iti], [0 max(veloc(2:end))*1.2],'b-', 'LineWidth',1)
        hold off
        
        
    end

end

