function hist2batch(eptrials, cells, bins)

% Heatmaps for all cells

%samplingrate
smplrt=length(eptrials(eptrials(:,4)==1,1))/max(eptrials(:,1));

%subset of trials
%eptrials = eptrials(eptrials(:,5)>20, :);


%FIRST DEAL WITH TIME (COMMON TO ALL CELLS). For comments, see hist2

%all time samples in two vectors
xt = eptrials(eptrials(:,4)==1 & eptrials(:,9)==1 & ismember(eptrials(:,10),1:7), 2);
yt = eptrials(eptrials(:,4)==1 & eptrials(:,9)==1 & ismember(eptrials(:,10),1:7), 3);

%evenly spaced bins of x and y coordinate ranges (incl pos - not just event -
%data)
xedges = linspace(min(eptrials(eptrials(:,9)==1 & ismember(eptrials(:,10),1:7),2)), max(eptrials(eptrials(:,9)==1 & ismember(eptrials(:,10),1:7),2)), bins);
yedges = linspace(min(eptrials(eptrials(:,9)==1 & ismember(eptrials(:,10),1:7),3)), max(eptrials(eptrials(:,9)==1 & ismember(eptrials(:,10),1:7),3)), bins);

[xnt, xbint] = histc(xt,xedges);
[ynt, ybint] = histc(yt,yedges);


        %THIS SECTION REMOVES PIXLES THAT WERE ONLY VISITED ON ONE TRIAL. IT IS NOT
        %APPROPRIATE FOR NON-TRIAL BASED DATA (CAN BE COMMENTED OUT).

        %matrix y-cords x-cords and trial number
        pixles_and_trials = [ybint xbint eptrials(eptrials(:,4)==1 & eptrials(:,9)==1 & ismember(eptrials(:,10),1:7), 6)]; 

        %get indices of unique rows (each pixle visited on each trial)
        [~, uni_indi, ~] = unique(pixles_and_trials, 'rows');
        pixle_per_trial = pixles_and_trials(uni_indi, 1:2);

        %get repeated rows of pixles_and_trials(indices, 1:2) (pixles that were 
        %visited on multiple trials)
        [~,~,n] = unique(pixle_per_trial, 'rows'); %see help unique
        hist_counts = hist(n, .5:1:(max(n)-.5)); %how many trials each pixle was visited on
        multi_trial_pixles = pixle_per_trial(ismember(n, find(hist_counts>1)), :); %pixles visited on >1 trials

        %index pixles_and_trials for pixles visited on multiple trials
        multivisit_pixle_index_space = logical(ismember([ybint xbint], multi_trial_pixles, 'rows'));

        %re-define set of visited x and y matrix coords to only include
        %pixles visited on multiple trials
        xbint = xbint(multivisit_pixle_index_space);
        ybint = ybint(multivisit_pixle_index_space);
                
        

xnbins = length(xedges);
ynbins = length(yedges);
xnbint = length(xedges);
ynbint = length(yedges);

if xnbint >= ynbint
    xyt = ybint*(xnbint) + xbint;
    indexshiftt = xnbint;
else
    xyt = xbint*(ynbint) + ybint;
    indexshiftt = ynbint;
end

xyunit = unique(xyt);
hstrest = histc(xyt,xyunit);

histmatt = zeros(xnbint, ynbint);

histmatt(xyunit-indexshiftt) = hstrest;
histmatt = histmatt'./smplrt;

%PREALLOCATE HERE
hold_matrices = nan(length(histmatt(:,1)), length(histmatt(1,:)), length(cells));


%LOOP THROUGH CELLS TO DEAL WITH SPIKES

for c = 1:length(cells)
    
    %all spike events in two vectors
    xs = eptrials(eptrials(:,4)==c & eptrials(:,9)==1 & ismember(eptrials(:,10),1:7), 2);
    ys = eptrials(eptrials(:,4)==c & eptrials(:,9)==1 & ismember(eptrials(:,10),1:7), 3);

    %filling xbin and ybin with firing rates. Last row is always 0, so we
    %remove it.

    %spikes
    [xns, xbins] = histc(xs,xedges);
    [yns, ybins] = histc(ys,yedges);
    
    %remove singleton pixles
    multivisit_pixle_index_event = logical(ismember([ybins xbins], multi_trial_pixles, 'rows'));
    xbins = xbins(multivisit_pixle_index_event);
    ybins = ybins(multivisit_pixle_index_event);


    %xbin, ybin zero for out of range values (see the help of histc) force this 
    %event to the first bins

    xbins(find(xbins == 0)) = 1;
    ybins(find(ybins == 0)) = 1;
    xbint(find(xbint == 0)) = 1;
    ybint(find(ybint == 0)) = 1;
    
    if xnbins >= ynbins
        xys = ybins*(xnbins) + xbins;
        indexshifts = xnbins;
    else
        xys = xbins*(ynbins) + ybins;
        indexshifts = ynbins;
    end

    xyunis = unique(xys);
    hstress = histc(xys,xyunis);
    
    histmats = zeros(xnbins, ynbins);

    histmats(xyunis-indexshifts) = hstress;
    histmats = histmats';
    

    %spikes
    histmats(1,1)=0;
    %time
    histmatt(1,1)=0;

    matrix = histmats./histmatt;

    if sum(isinf(matrix(:)))>0
        %disp('Warning: firing rate bins containing infinite values were set as NaN')
        matrix(isinf(matrix))=NaN;
    end

    size_fix = matrix;
    size_fix(~isnan(size_fix))=1;
    
    %convolve
    mask = [1 3 1; 3 4 3; 1 3 1]./20;
    
    matrix = conv2nan(matrix, mask, 'same');
    matrix = conv2nan(matrix, mask, 'same');
    matrix = conv2nan(matrix, mask, 'same');
    matrix = conv2nan(matrix, mask, 'same');
    matrix = conv2nan(matrix, mask, 'same');
    
    %remove convolve bloat
    matrix = matrix.*size_fix;
    
    %plot
    figure; 
    pcolor(xedges,yedges,matrix); 
    colorbar; 
    axis square tight;
    set(gca,'xdir','reverse')
    shading flat
    caxis([0 max(matrix(:))*.7])
    title(['Cell ',num2str(cells(c))],'fontsize', 16)
    
    %save figures
    saveas(gcf,fullfile('/Users/ampm/Desktop/temp',[num2str(cells(c)),'.fig']),'fig')
    saveas(gcf,fullfile('/Users/ampm/Desktop/temp',[num2str(cells(c)),'.jpg']),'jpg')
    
end

end