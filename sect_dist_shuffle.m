function [sdssa mean_sdssa lowrng_sdssa highrng_sdssa sds1, sds2, sds3] = sect_dist_shuffle(num_shuffles, stages, varargin)
%calculates the sec distances after each of num_shuffles of each of the 
%stages in vargargin
tic
sdssa = nan(num_shuffles,5,3);

%for each stage variable
for stage = stages
    stage
    stg_mtx = varargin{stage};
    stg_cols = [];
    sds = nan(size(stg_mtx,3), 5);
    clust_count = 0;
    
    %reshape variable from 3d to 2d
    for page = 1:size(stg_mtx,3)
        stg_cols = [stg_cols; stg_mtx(:,:,page)];
    end
    
    shuffle_count = 0;
    for shuffle = 1:num_shuffles
        shuffle_count = shuffle_count+1;
    
        %randomize LR trials within a session (index first cluster in sesh)
        for sesh = 1:max(stg_cols(:,8))
            perm = randperm(sum(stg_cols(:,8)==sesh & stg_cols(:,7)==1 & ~isnan(stg_cols(:,6)) & stg_cols(:,9)==1));
            orig = stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==1 & ~isnan(stg_cols(:,6)) & stg_cols(:,9)==1, 6);
            shuf = orig(perm);

            %shuffled L's R's and original NaN's
            shuf_idx = nan(size(stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==1  & stg_cols(:,9)==1, 6)));  
            shuf_idx(~isnan(stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==1 & stg_cols(:,9)==1, 6))) = shuf;

            %non_shuffle
            %shuf_idx(~isnan(stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==1 & stg_cols(:,9)==1,6))) = orig;
            %test = [stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==1 & ~isnan(stg_cols(:,6)) & stg_cols(:,9)==1, 6) shuf_idx];

            %use shuffled LR trial order to calculate sec dists for each
            %cluster in session
            for cluster = 1:max(stg_cols(stg_cols(:,8)==sesh, 7))

                clust_count = clust_count+1;
                sect_rates = stg_cols(stg_cols(:,8)==sesh & stg_cols(:,7)==cluster & stg_cols(:,9)==1, 1:5);

                %index for left and right clust/sector/session rates
                Ls = sect_rates(shuf_idx==1,:);
                Rs = sect_rates(shuf_idx==2,:);

                %section means
                lmeans = nanmean(Ls);
                rmeans = nanmean(Rs);
                
                %{
                %sections stds
                lstds = nanstd(Ls);
                rstds = nanstd(Rs);
                
                %pooled std
                pstd_app = nan(1,5);
                for i = 1:5
                    %[~,~,pstd_app_temp] = pooledmeanstd(length(Ls(~isnan(Ls(:,i)),i)),lmeans(i),lstds(i),length(Rs(~isnan(Rs(:,i)),i)),rmeans(i),rstds(i));
                    [~,~,~, pstd] = ttest2(Ls(~isnan(Ls(:,i)),i),Rs(~isnan(Rs(:,i)),i));
                    pstd_app_temp = (pstd.sd);
                    
                    pstd_app(i) = pstd_app_temp;
                end
                %}
                %load standard difference
                sds(clust_count, :) = abs(lmeans-rmeans);%./pstd_app;

            end
        end
        
        switch stage
            case 1
                sds1 = sds;
            case 2
                sds2 = sds;
            case 3
                sds3 = sds;
        end
        
        %load shuffle iteration mean
        sdssa(shuffle, :, stage) = nanmean(sds);
        
    end
    
    
end

%calculate mean and 99% range
mean_sdssa = mean(sdssa);
sort_sdssa = sort(sdssa);
if floor(size(sort_sdssa,1)*.005) >0
    lowrng_sdssa = sort_sdssa(floor(size(sort_sdssa,1)*.005), :, :);
else 
    lowrng_sdssa = sort_sdssa(1, :, :);
end
highrng_sdssa = sort_sdssa(ceil(size(sort_sdssa,1)*.995), :, :);
toc
end