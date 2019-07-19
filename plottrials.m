function plottrials(eptrials, section_boundaries, varargin)
%function plottrials(eptrials, cell_range)
%plot the trajectory for all trials, or just a specified range


%check inputs
switch nargin
    case 0
        error(message('too few inputs'))
    case 1
        error(message('too few inputs'))
    case 2
    case 3
        trials = varargin{1};           
    otherwise
        %error(message('too many inputs'))
end


%what trials to plot
if exist('trials', 'var')
    range = trials;
else
    range = 1:1:(max(eptrials(:,6)));
end

for trial = range 
    
    local_eptrials = eptrials(eptrials(:,6)==trial,:);
    trl_min_time = min(local_eptrials(:,1));
    number_timestamps_trl = length(local_eptrials(:,1));
    time_correction_vector = ones(number_timestamps_trl,1)*trl_min_time;
    local_eptrials(:,1) = (local_eptrials(:,1) - time_correction_vector);
    
    figure
    hold on
    
    %plot entire trajectory in black
    plot3(local_eptrials(:,2), local_eptrials(:,3), local_eptrials(:,1), 'k')
    
    %plot portion of trajectory where light is on in red
    plot3(local_eptrials(ismember(local_eptrials(:,13), [1 2]),2), local_eptrials(ismember(local_eptrials(:,13), [1 2]),3), local_eptrials(ismember(local_eptrials(:,13), [1 2]),1), 'r.')%, 'LineWidth', 2)
    
    %plot portion of trajectory where fake light is "on" in cyan
    plot3(local_eptrials(ismember(local_eptrials(:,13), 3),2), local_eptrials(ismember(local_eptrials(:,13), 3),3), local_eptrials(ismember(local_eptrials(:,13), 3),1), 'c.')%, 'LineWidth', 2)
    
    %plot portion of trajectory categorized as ITI in magenta
    plot3(local_eptrials(local_eptrials(:,10)==8,2), local_eptrials(local_eptrials(:,10)==8,3), local_eptrials(local_eptrials(:,10)==8,1), 'm.')%, 'LineWidth', 2)
     
    %plot portion of trajectory categorized as carry in green
    plot3(local_eptrials(local_eptrials(:,10)==9,2), local_eptrials(local_eptrials(:,10)==9,3), local_eptrials(local_eptrials(:,10)==9,1), 'g.')%, 'LineWidth', 2)
    
    %blue dot at reward flag
    plot3(local_eptrials(ismember(local_eptrials(:,5), [1 2]),2), local_eptrials(ismember(local_eptrials(:,5), [1 2]),3), local_eptrials(ismember(local_eptrials(:,5), [1 2]),1), 'b.', 'markersize', 20)
    
    
    
    %maze boundaries
    sections(eptrials, section_boundaries)
    
    %plot start and stop icons
    try
        plot3(local_eptrials(local_eptrials(:,1)==min(local_eptrials(:,1)), 2), local_eptrials(local_eptrials(:,1)==min(local_eptrials(:,1)), 3), min(local_eptrials(:,1)), 'g*', 'markersize', 20)
        plot3(local_eptrials(local_eptrials(:,1)==max(local_eptrials(:,1)), 2), local_eptrials(local_eptrials(:,1)==max(local_eptrials(:,1)), 3), max(local_eptrials(:,1)), 'r*', 'markersize', 20)
    catch
        trial
    end
    %title
    
    if sum(isnan(local_eptrials(:,7)))>sum(local_eptrials(:,7)==mode(local_eptrials(:,7)))
        
        %left or right trial type (for plots)
        if mode(local_eptrials(:, 8))==1
            type = 'Left';
        elseif mode(local_eptrials(:, 8))==2
            type = 'Right';
        else
            type = 'Unknown_L/R';
        end
    
        title(['Trial ',num2str(trial), ' ',num2str(type)],'fontsize', 16)
    
    else
    
        %left or right trial type (for plots)
        if mode(local_eptrials(:, 7))==1
            type = 'Left';
        elseif mode(local_eptrials(:, 7))==2
            type = 'Right';
        else
            type = 'Unknown_L/R';
        end
    
        %correct or error trial type (for plots)
        if mode(local_eptrials(:, 9))==1
            accuracy = 'Correct';
        elseif mode(local_eptrials(:, 9))==0
            accuracy = 'Error';
        else
            accuracy = 'UnknownAccuracy';
        end
        
        title(['Trial ',num2str(trial), ' ',num2str(type), ' ',num2str(accuracy)],'fontsize', 16)
    
    end
    
    %reverse z (time) plot so that time increase as rat moves away from
    %viewer default angle
    set(gca,'zdir','reverse')

end