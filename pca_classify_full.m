function [class_pca_prop, class_pca, shuf_out, comb_safters, ID_idx] = pca_classify_full(stage_numbers, LR, evt, shuffs)
%function calculates firing rates at multiple times on the maze and
%then attempts to classify the times based on the firing rates
%
%times on maze determined by ALL_keytimes

%calculate firing rates
[safter_L, IDs_L] = ALL_keytimes_LR_full(stage_numbers, 1, evt);
[safter_R, IDs_R] = ALL_keytimes_LR_full(stage_numbers, 2, evt);

if LR == 1
    high_oneside = max(IDs_L);
    IDs_R = IDs_R + repmat(high_oneside, size(IDs_R));
end

ID_idx = [IDs_L; IDs_R];

%calculate pca
comb_safters = [safter_L; safter_R];

%standardize
comb_safters = comb_safters - repmat(mean(comb_safters), size(comb_safters,1), 1);
stds = std(comb_safters);
stds(stds==0) = 1;
comb_safters = comb_safters./repmat(stds, size(comb_safters,1), 1);


%comb_safters = comb_safters(randperm(size(comb_safters,1)), :);
[~, pc_vectors, ~, ~, variance_explained] = pca(comb_safters);



%classify based on a subset of dimensions (dims)
dims = 1:100;

%classify either all data (classify one trial at a time), or train first half test seconed half
class_pca = nan(size(pc_vectors,1),1);
for i = 1:size(pc_vectors,1)
    class_pca(i) = classify(pc_vectors(i, dims), pc_vectors(setdiff(1:size(pc_vectors,1),i), dims), ID_idx(setdiff(1:size(pc_vectors,1),i)));
end
   
%calculate success
class_pca_prop = sum(class_pca == ID_idx)/length(class_pca);

%FIGURES (hardish coded)
%

%PCA 2d plots
%
components = [1 2]; 
figure; hold on
colors = [255 225 102; 255 178 102; 216.7500 82.8750 24.9900; 142 0 0; 102 205 255; 51 153 255;  0  113.9850  188.9550; 50 0 150];

if numel(components) == 2
    for id = unique(ID_idx)'
        plot(pc_vectors(ID_idx==id, components(1)), pc_vectors(ID_idx==id, components(2)), '.', 'Color', colors(id,:)./255, 'MarkerSize', 30);
    end
    borderx = (max(pc_vectors(:,components(1))) - min(pc_vectors(:,components(1))))/10;
    bordery = (max(pc_vectors(:,components(2))) - min(pc_vectors(:,components(2))))/10;
    axis([min(pc_vectors(:,components(1)))-borderx max(pc_vectors(:,components(1)))+borderx min(pc_vectors(:,components(2)))-bordery max(pc_vectors(:,components(2)))+bordery])
elseif numel(components) == 3
    for id = unique(ID_idx)'
        plot3(pc_vectors(ID_idx==id, components(1)), pc_vectors(ID_idx==id, components(2)), pc_vectors(ID_idx==id, components(3)), '.', 'Color', colors(id,:)./255, 'MarkerSize', 30);
    end    
end
box off; set(gca,'TickLength',[0, 0]);
hold off
%}

%variance explained plot
%
%
figure;
var_exp = cumsum(variance_explained);
plot(var_exp)
box off; set(gca,'TickLength',[0, 0]);
%variance explained by dims of pca
var_exp(max(dims))

%}

%classification success bar plot
%
figure; hold on

%event-specific classification success rates
bar_in = nan(1, length(unique(ID_idx)));
for id = unique(ID_idx)'
    bar_in(id) = sum(class_pca(ID_idx==id) == ID_idx(ID_idx==id));
end
%barchart

if LR == 1
    bar_in = mean(reshape(bar_in, length(bar_in)/2, 2),2);
end

bar(bar_in./length(class_pca(ID_idx==1)))

ylim([0 1.05])
box off; set(gca,'TickLength',[0, 0]);


%shuffle
%
if LR == 0
    shuf_out = nan(shuffs, length(unique(ID_idx))+1);
elseif LR == 1
    shuf_out = nan(shuffs, (length(unique(ID_idx))/2)+1);
end

warning('off','all') %pca cries about dimensionality during shuffles. 
                     %That's expected. The shuffle is creating a flat 
                     %(i.e. random) dimension.

for shuf = 1:shuffs

    %shuffle index matrix
    shufled_comb_safters = nan(size(comb_safters)); 
    
    %if LR == 0 %shuffle each cell's timebin rates between all timebins
        for i = 1:size(comb_safters,2)
            shufled_comb_safters(:,i) = comb_safters(randperm(size(comb_safters,1)), i);
        end
    %elseif LR == 1 %shuffle each cell's timebin rates between the LR timebins of that event
        %for column = 1:size(shuf_idx,2)
            %shuf_idx(:,column) = 1:size(shuf_idx,1);
            %for event = 1:length(unique(train_IDs_L))
                %shuf_hold = shuf_idx(ismember(ID_idx, [event, event+high_oneside]), column);
                %shuf_idx(ismember(ID_idx, [event, event+high_oneside]), column) = shuf_hold(randperm(size(shuf_hold,1)));
            %end
        %end
    %end

    %pca shuffle
    [~, pc_vectors_shuf, ~, ~, ~] = pca(shufled_comb_safters);

    %classify and report success
    class_pca_shuf = nan(size(pc_vectors_shuf,1),1);
    for i = 1:size(pc_vectors_shuf,1)
        class_pca_shuf(i) = classify(pc_vectors_shuf(i, dims), pc_vectors_shuf(setdiff(1:size(pc_vectors_shuf,1),i), dims), ID_idx(setdiff(1:size(pc_vectors_shuf,1),i)));
    end
    shuf_out(shuf, 1) = sum(class_pca_shuf == ID_idx)/length(class_pca_shuf);
    
    
    %report event-specific success
    if LR == 0
        for col = 2:size(shuf_out,2)
            shuf_out(shuf, col) = sum(class_pca_shuf(ID_idx==col-1) == ID_idx(ID_idx==col-1))./sum(ID_idx==1);
        end
    elseif LR == 1
        for col = 2:size(shuf_out,2)
            %train_idx = ismember(train_IDs,[col col+high_oneside]); 
            %test_idx = ismember(test_IDs,[col col+high_oneside]);
            shuf_out(shuf, col) = mean([sum(class_pca_shuf(ID_idx==col-1) == ID_idx(ID_idx==col-1))./sum(ID_idx==1), sum(class_pca_shuf(ID_idx==col-1+high_oneside) == ID_idx(ID_idx==col-1+high_oneside))./sum(ID_idx==1)]);
        end
    end
end

warning('on','all')

%overlay shuffle output on barchart
shuf_out_sorted = sort(shuf_out);

if length(mean(shuf_out_sorted(:,2:end))) == 1 
    plot([.25 size(shuf_out,2)-.25], [mean(shuf_out_sorted(:,2:end)) mean(shuf_out_sorted(:,2:end))], 'k-') %mean shuffles
    plot([.25 size(shuf_out,2)-.25], [shuf_out_sorted(floor(size(shuf_out_sorted,1)*.025), 2:end) shuf_out_sorted(floor(size(shuf_out_sorted,1)*.025), 2:end)], 'k--') %low bound
    plot([.25 size(shuf_out,2)-.25], [shuf_out_sorted(ceil(size(shuf_out_sorted,1)*.975), 2:end) shuf_out_sorted(ceil(size(shuf_out_sorted,1)*.975), 2:end)], 'k--') %high bound
else    
    plot([.25 1:1:size(shuf_out,2)-1 size(shuf_out,2)-.25], [nan mean(shuf_out_sorted(:,2:end)) nan], 'k-') %mean shuffles
    plot([.25 1:1:size(shuf_out,2)-1  size(shuf_out,2)-.25], [nan shuf_out_sorted(floor(size(shuf_out_sorted,1)*.025), 2:end) nan], 'k--') %low bound
    plot([.25 1:1:size(shuf_out,2)-1  size(shuf_out,2)-.25], [nan shuf_out_sorted(ceil(size(shuf_out_sorted,1)*.975), 2:end) nan], 'k--') %high bound
end

end