% Script used to obtain the results for the clustering project manuscript
%%%

%%author__ = Shanna Kulik & Jolanda Derks
%%contact__ = s.kulik@amsterdamumc.nl
%%date__ = 2020
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2021/12/2

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% David Ferreira (2021). Quantiles (https://www.mathworks.com/matlabcentral/fileexchange/70279-quantiles), MATLAB Central File Exchange. Retrieved December 9, 2021.

%%% Toolboxes 
% Brain Connectivity Toolbox https://sites.google.com/site/bctnet/ 
% version 2017-15-01



%% 1. Global clustering differences between HC & pt
load('/path/to/AECc_load_matrices.mat')

norm_cl_HC = mean(cl_HC,4)./mean(shuf_rand_matrix_cl_HC,4);
norm_cl_pt = mean(cl_NO,4)./mean(shuf_rand_matrix_cl_NO,4);

mean_norm_cl_HC = squeeze(mean(mean(norm_cl_HC,3),2));
mean_norm_cl_NO = squeeze(mean(mean(norm_cl_pt,3),2));

% Test global clustering HC vs pt
[p, h, stats] = ranksum(mean_norm_cl_HC,mean_norm_cl_NO);                       % p=0.0384, U=2902 --> new 2020: p>0.01, U=2353

[Q,IQR] = quartile(mean_norm_cl_NO);
med = median(mean_norm_cl_NO);
[Q,IQR] = quartile(mean_norm_cl_HC);
med = median(mean_norm_cl_HC);

figure()
boxplot([mean_norm_cl_NO' mean_norm_cl_HC'],[zeros(1,71) ones(1,53)])
xlabel(['Mean glob clus NO', 'Mean glob clus HC'])
ylabel('Normalized (shuffled) global clustering value')

%% 2. Tumor vs non-tumor pt

subj_pt = 1:71;
epochs_pt = 1:15; 
sub_epochs = 1:4;

for sub = 1:numel(subj_pt)
    for ep = 1:numel(epochs_pt)   
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(AEC_theta_NO_zero(sub,ep,sub_ep,:,:));
            cl_pt_theta(sub,ep,sub_ep,:) = clustering_coef_wu(matrix);              % size matrix = 71x15x4x78 (clustering per region, per epoch, per patient)
        end
    end
end

% Norm-scores AEC & clustering patients
    
for sub = 1:numel(subj_pt)
    for ep = 1:numel(epochs_pt)
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(AEC_theta_NO_zero(sub,ep,sub_ep,:,:));
            matrix_norm = weight_conversion(matrix, 'normalize');
            AEC_theta_norm_pt(sub,ep,sub_ep,:,:) = matrix_norm;                             % Normalized PLI matrix of each patient, size=71x15x4x78x78
            matrix_norm_AEC = squeeze(AEC_theta_norm_pt(sub,ep,sub_ep,:,:));
            clus_theta_norm_pt(sub,ep,sub_ep,:) = clustering_coef_wu(matrix_norm_AEC);        % Normalized local clus matrix of each patient, size=71x15x4x78
        end
    end
end

av_roi_cl_norm_pt = squeeze(mean(mean(clus_theta_norm_pt,3),2));                              % Mean normalized local clus per patient, averaged over epochs, size 71x78 

% tumor roi and  non-tumor roi clustering (weighted for each tumor roi)

% import data nr of voxels per AAL-MNI region
cd('/path/to/masks')                                

files = dir('*_AAL.txt');                                                                 % 71 files 
subj = 1:numel(files);

tumor_clus = zeros(numel(subj),1);

for k = 1:numel(subj)  
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL = importdata(fullFileName);
    tumor_clus_pt(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*av_roi_cl_norm_pt(k,:)');     % Clustering of weighted tumor regions taken together per patient, size 71x1
    non_tumor_voxels = find(voxels_AAL(1:78,1)==0);
    non_tumor_clus_pt(k,:) = mean(av_roi_cl_norm_pt(k,non_tumor_voxels));                              % Clustering of average of all non tumor regions per patient, size 71x1
end

% Test difference between tumor and non tumor clustering in patients 

[p,h,stats]=signrank(tumor_clus_pt, non_tumor_clus_pt) % p = 0.2814 W=1466 --> new 2020: p=0.2865, W=1464

[Q,IQR] = quartile(tumor_clus_pt);
med = median(tumor_clus_pt);
[Q,IQR] = quartile(non_tumor_clus_pt);
med = median(non_tumor_clus_pt);

figure()
boxplot([tumor_clus_pt non_tumor_clus_pt],[zeros(1,71) ones(1,71)])
xlabel(['tumor clus', 'non-tumor clus'])
ylabel('Local normalized clustering value')
%% 3. Z-score tumor vs non-tumor
subj_HC = 1:53;
epochs_HC = 1:13; 
    
% Norm scores HC

for sub = 1:numel(subj_HC)
    for ep = 1:numel(epochs_HC)   
        for sub_ep = 1:numel(sub_epochs)
        matrix = squeeze(AEC_theta_HC_zero2(sub,ep,sub_ep,:,:));
        matrix_norm = weight_conversion(matrix, 'normalize');
        AEC_theta_norm_HC(sub,ep,sub_ep,:,:) = matrix_norm;                     % Normalized PLI matrix of each HC, size=71x15x4x78x78
        matrix2 = squeeze(AEC_theta_norm_HC(sub,ep,sub_ep,:,:));
        clus_theta_norm_HC(sub,ep,sub_ep,:) = clustering_coef_wu(matrix2);      % Normalized local clus matrix of each HC, size=65x15x4x78
        end
    end
end

av_cl_norm_HC = squeeze(mean(mean(mean(clus_theta_norm_HC,4),3),2));            % Mean normalized local clus per HC, size=65x1
av_roi_cl_norm_HC = squeeze(mean(mean(mean(clus_theta_norm_HC))));              % Mean normalized local clus per roi, size=78x1
av_roi_cl_norm_HC2 = squeeze(mean(mean(mean(clus_theta_norm_HC,3),2),1));

% HC tumor roi and  non-tumor roi clustering (weighted for each tumor roi)

cd('/path/to/masks')

tumor_clus_HC = zeros(numel(subj),1);

for k = 1:numel(subj)  
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL = importdata(fullFileName);
    tumor_clus_HC(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*av_roi_cl_norm_HC);    % Clustering of weighted tumor regions taken together per HC, size 71x1
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    non_tumor_clus_HC(k,:) = mean(av_roi_cl_norm_HC(non_tumor_voxels_AAL));                     % Clustering of average of all non tumor regions per HC, size 71x1
end

% z-scores obv HC clustering 

av_cl_roi_HC=squeeze(mean(mean(clus_theta_norm_HC,3),2));                     % Local clustering for each region of each HC - 65x78
sd_av_cl_roi_HC = std(av_cl_roi_HC,0,1)';                                   % 78x1

for i = 1:71
    z_score(:,i) = (av_roi_cl_norm_pt(i,:)' - av_roi_cl_norm_HC)./sd_av_cl_roi_HC;  % clustering of each patient - clustering of HC regions/sd of hc regions
end 

z_score = z_score';

for k = 1:numel(subj) 
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL = importdata(fullFileName);
    z_cl_roi_tumor(k,:) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*z_score(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    z_cl_roi_non_tumor(k,:) = mean(z_score(k,non_tumor_voxels_AAL));
end

% Test difference between tumor and non-tumor regions with z-scores
[p h stats] = signrank(z_cl_roi_tumor,z_cl_roi_non_tumor); % 

[Q,IQR] = quartile(z_cl_roi_tumor);
med = median(z_cl_roi_tumor);
[Q,IQR] = quartile(z_cl_roi_non_tumor);
med = median(z_cl_roi_non_tumor);

figure()
boxplot([z_cl_roi_tumor,z_cl_roi_non_tumor],[zeros(1,71) ones(1,71)])
xlabel(['intrinsic tumor clustering', 'intrinsic non-tumor clustering'])
ylabel('local normalized clustering value')

%% 4. Euclidian distance whole group

load('/path/to/AECc_euclidian.mat')

res_euc_dis = reshape(euc_distance, 5538, [ ]);
res_hyp_dis = reshape(hyperbolic_distance,5538, [ ]);

z_score_flip = z_score';
res_cl_pt_roi_z = reshape(z_score_flip,5538,[ ]);

res_euc_dis_z(res_euc_dis<=0) = NaN;
res_hyp_dis_z(res_hyp_dis<=0) = NaN;

res_euc_dis(res_euc_dis<=0) = NaN;
res_hyp_dis(res_hyp_dis<=0) = NaN;

[rho_euc_z, pval_euc_z] = corr(res_euc_dis, res_cl_pt_roi_z,'rows','complete','type', 'Spearman'); % rho=0.0230, p=0.0866 --> new 2020: rho=0.0184, p=0.1718
[rho_hyp_z, pval_hyp_z] = corr(res_hyp_dis, res_cl_pt_roi_z,'rows','complete','type', 'Spearman'); % rho=-.0088, p=0.5128

figure(); 
subplot(2,1,1)
scatter(res_euc_dis, res_cl_pt_roi_z, 'r')
title('Euclidian distance in relation to local clustering')
lsline
subplot(2,1,2)
scatter(res_hyp_dis, res_cl_pt_roi_z, 'k')
title('Hyperbolic distance in relation to local clustering')
lsline

%% Linear mixed model

grouping = ones(71*78,1);

for i = 1:71
    grouping((i-1)*78+1:i*78)= i;
end
LME = fitlme(table(res_cl_pt_roi_z,res_euc_dis,grouping),'res_euc_dis~res_cl_pt_roi_z+(1|grouping)+(res_cl_pt_roi_z-1|grouping)');

figure()
scatter(res_euc_dis, res_cl_pt_roi_z, 'k'); hold on
scatter(res_euc_dis(1:78), res_cl_pt_roi_z(1:78), 'r', 'filled')
scatter(res_euc_dis(2*78-1:3*78), res_cl_pt_roi_z(2*78-1:3*78), 'g', 'filled')
scatter(res_euc_dis(5*78-1:6*78), res_cl_pt_roi_z(5*78-1:6*78), 'b', 'filled')
%% 5. High vs low tumor z-scores   

mean_norm_glob_clus_pt = squeeze(mean(mean(norm_cl_pt,3),2));                % mean global clus per patient, size=71x1
norm_clus_pt_high = find(mean_norm_glob_clus_pt > 1.001436);                 % n=18  old: 1.001386
norm_clus_pt_low = find(mean_norm_glob_clus_pt < 1.001436);                  % n=53

[p h stats] = ranksum(z_cl_roi_tumor(norm_clus_pt_low),z_cl_roi_tumor(norm_clus_pt_high)); % 2-3-2020 new: p<0.001, U=1525

figure()
boxplot([z_cl_roi_tumor(norm_clus_pt_low)' z_cl_roi_tumor(norm_clus_pt_high)'],[zeros(1,53) ones(1,18)])
xlabel(['tumor clustering low group', 'Tumor clustering high group'])
ylabel('local clustering value')

%% 6. High vs low non-tumor z-scores  

[p h stats] = ranksum(z_cl_roi_non_tumor(norm_clus_pt_low), z_cl_roi_non_tumor(norm_clus_pt_high)); %2-3-2020 new: p<0.001, U=1508

[Q,IQR] = quartile(z_cl_roi_non_tumor(norm_clus_pt_low));
med = median(z_cl_roi_non_tumor(norm_clus_pt_low));
[Q,IQR] = quartile(z_cl_roi_non_tumor(norm_clus_pt_high));
med = median(z_cl_roi_non_tumor(norm_clus_pt_high));

figure()
boxplot([z_cl_roi_non_tumor(norm_clus_pt_low)' z_cl_roi_non_tumor(norm_clus_pt_high)'],[zeros(1,53) ones(1,18)])
xlabel(['non-tumor clustering low group', 'non-tumor clustering high group'])
ylabel('local clustering value')

%% 7. Euclidian distance high & low

high_group = [20 70 54 64 7 40 12 34 35 19 21 3 6 26 4 37 27 23]; % matrix positions for high clustering group
euc_dis_low = euc_distance;
hyp_dis_low = hyperbolic_distance;
euc_dis_low(high_group,:)=[];
hyp_dis_low(high_group,:)=[];

euc_dis_high = euc_distance(high_group,:);
hyp_dis_high = hyperbolic_distance(high_group,:);

res_euc_dis_tumor_high = reshape(euc_dis_high, 1404,[ ]);
res_euc_dis_tumor_low = reshape(euc_dis_low, 4134,[ ]);
res_hyp_dis_tumor_high = reshape(hyp_dis_high, 1404,[ ]);
res_hyp_dis_tumor_low = reshape(hyp_dis_low, 4134,[ ]);

z_score_high = z_score_flip(:,(norm_clus_pt_high));
z_score_low = z_score_flip(:,(norm_clus_pt_low));

res_z_score_high = reshape(z_score_high, 1404, [ ]);
res_z_score_low = reshape(z_score_low, 4134, [ ]);

[rho_euc_z_high, pval_euc_z_high] = corr(res_euc_dis_tumor_high, res_z_score_high, 'type', 'Spearman'); % rho= 0.0145, p=0.5871
[rho_euc_z_low, pval_euc_z_low] = corr(res_euc_dis_tumor_low, res_z_score_low, 'type', 'Spearman'); % rho=-0.0098 , p=0.5292
[rho_hyp_z_high, pval_hyp_z_high] = corr(res_hyp_dis_tumor_high, res_z_score_high, 'type', 'Spearman'); % rho= -0.0162 , p= 0.5445
[rho_hyp_z_low, pval_hyp_z_low] = corr(res_hyp_dis_tumor_low, res_z_score_low, 'type', 'Spearman'); % rho= -0.0195, p=0.2096

figure(); 
subplot(2,2,1)
scatter(res_euc_dis_tumor_high,res_z_score_high, 'r')
title('high clustering group - euclidian distance')
lsline
subplot(2,2,2)
scatter(res_euc_dis_tumor_low,res_z_score_low, 'k')
title('low clustering group - euclidian distance')
lsline
subplot(2,2,3)
scatter(res_hyp_dis_tumor_high,res_z_score_high, 'r')
title('high clustering group - hyperbolic distance')
lsline
subplot(2,2,4)
scatter(res_hyp_dis_tumor_low,res_z_score_low, 'k')
title('low clustering group - hyperbolic distance')
lsline

%% Linear mixed model for high and low groups

grouping_high = ones(18*78,1);
grouping_low = ones(53*78,1);

for i = 1:18
    grouping_high((i-1)*78+1:i*78)= i;
end

for i = 1:53
    grouping_low((i-1)*78+1:i*78)= i;
end

LME_high = fitlme(table(res_z_score_high,res_euc_dis_tumor_high,grouping_high),'res_euc_dis_tumor_high~res_z_score_high+(1|grouping_high)+(res_z_score_high-1|grouping_high)');
LME_low = fitlme(table(res_z_score_low,res_euc_dis_tumor_low,grouping_low),'res_euc_dis_tumor_low~res_z_score_low+(1|grouping_low)+(res_z_score_low-1|grouping_low)');

figure()
scatter(res_euc_dis_tumor_high, res_z_score_high, 'k'); hold on
scatter(res_euc_dis_tumor_high(1:78), res_z_score_high(1:78), 'r', 'filled')
scatter(res_euc_dis_tumor_high(2*78-1:3*78), res_z_score_high(2*78-1:3*78), 'g', 'filled')
scatter(res_euc_dis_tumor_high(5*78-1:6*78), res_z_score_high(5*78-1:6*78), 'b', 'filled')

figure()
scatter(res_euc_dis_tumor_low, res_z_score_low, 'k'); hold on
scatter(res_euc_dis_tumor_low(1:78), res_z_score_low(1:78), 'r', 'filled')
scatter(res_euc_dis_tumor_low(2*78-1:3*78), res_z_score_low(2*78-1:3*78), 'g', 'filled')
scatter(res_euc_dis_tumor_low(5*78-1:6*78), res_z_score_low(5*78-1:6*78), 'b', 'filled')

%% Individual relations euclidian distance with global clustering

for i=1:71
    res_euc_dis_ind(:,i) = reshape(euc_distance(i,:),78, [ ]);
    res_hyp_dis_ind(:,i) = reshape(hyperbolic_distance(i,:),78, [ ]);
    res_cl_pt_roi_z_ind(:,i) = reshape(z_score_flip(:,i),78,[ ]);
    
    res_euc_dis_z_ind(res_euc_dis_ind<=0) = NaN;
    res_hyp_dis_z_ind(res_hyp_dis_ind<=0) = NaN;
    
    res_euc_dis_ind(res_euc_dis_ind<=0) = NaN;
    res_hyp_dis_ind(res_hyp_dis_ind<=0) = NaN;
    
    [rho_euc_z_ind(i), pval_euc_z_ind(i)] = corr(res_euc_dis_ind(:,i), res_cl_pt_roi_z_ind(:,i),'rows','complete','type', 'Spearman'); % rho=0.0230, p=0.0866
    [rho_hyp_z_ind(i), pval_hyp_z_ind(i)] = corr(res_hyp_dis_ind(:,i), res_cl_pt_roi_z_ind(:,i),'rows','complete','type', 'Spearman'); % rho=-.0088, p=0.5128
end

sig_euc_dis_ind=find(pval_euc_z_ind <0.05);

for i=1:20
    figure(); 
    scatter(res_euc_dis_ind(:,sig_euc_dis_ind(i)),res_cl_pt_roi_z_ind(:,sig_euc_dis_ind(i)))
    title('Individual euclidian distance in relation to local clustering')
    lsline
end

