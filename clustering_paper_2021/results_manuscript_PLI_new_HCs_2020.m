% Script used to obtain the results for the clustering project manuscript
%%%

%%author__ = Shanna Kulik & Jolanda Derks
%%contact__ = s.kulik@amsterdamumc.nl
%%date__ = 2020
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2021/12/27

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% David Ferreira (2021). Quantiles (https://www.mathworks.com/matlabcentral/fileexchange/70279-quantiles), MATLAB Central File Exchange. Retrieved December 9, 2021.
% https://github.com/RainCloudPlots/RainCloudPlots/blob/master/tutorial_matlab/raincloud_plot.m

%%% Toolboxes 
% Brain Connectivity Toolbox https://sites.google.com/site/bctnet/ 
% version 2017-15-01

%% Load & manipulate matrices
load('path/to/PLI_load_matrices.mat')

PLI_NO_theta = PLI_theta;
delete = [33 36 39 41 43 44 46 48 50 52 54 56 59 60 67 71 73 74 76 79];         % not included in this project
PLI_NO_theta(delete,:,:,:,:) = []; 

delete_HC = [9 27 28 35 37 38 44 45 53 55 60 64];
PLI_HC_theta(delete_HC,:,:,:,:) = [];

%% 1. Global clustering differences between HC & pt --> fig 1
cl_HC_all(11:12,:,:,:)=cl_HC2;
shuf_rand_matrix_cl_HC_all(11:12,:,:,:,:)=shuf_rand_matrix_cl_HC2;

norm_cl_pt = mean(cl_NO,4)./mean(shuf_rand_matrix_cl_NO,4);
norm_cl_NO = squeeze(mean(mean(norm_cl_pt,3),2));

norm_cl_hc = mean(cl_HC_all,4)./mean(shuf_rand_matrix_cl_HC_all,4);
norm_cl_HC = squeeze(mean(mean(norm_cl_hc,3),2));

norm_cl_HC(delete_HC,:) = []; 
cl_HC_all(delete_HC,:,:,:)=[];
shuf_rand_matrix_cl_HC_all(delete_HC,:,:,:)=[];

% Test global clustering HC vs pt
[p_clus_theta_new, tbl_clus_theta_new, stats_clus_theta_new] = ranksum(norm_cl_NO, norm_cl_HC); % p=0.0023 U=5042 --> nieuw: p=0.0024 U=5040
[Q,IQR] = quartile(norm_cl_NO);
med = median(norm_cl_NO);
[Q,IQR] = quartile(norm_cl_HC);
med = median(norm_cl_HC);

X=[1,1.5];
figure()
boxplot([norm_cl_NO' norm_cl_HC'],[zeros(1,71) ones(1,53)], 'positions', X)
xlabel(['Pt', 'HC'])
ylabel('Global normalized clustering value')

% Raincloud plot
figure; 
raincloud_plot(norm_cl_NO,'box_on', 1, 'dot_dodge_amounts', 0.1)

%% 2&3. Average theta FC & power

rel_pow_theta_HC_all(11:12,:,:,:)=rel_pow_theta_HC2;
mean_rel_pow_HC_theta_all(11:12,:,:,:)=mean_rel_pow_HC_theta2;
rel_pow_theta(delete,:,:,:) = []; 

mean_rel_pow_HC_theta_all(delete_HC,:)=[];
rel_pow_NO_theta = rel_pow_theta;
mean_rel_pow_pt_theta(delete,:) = [];

mean_PLI_HC = mean(mean(mean(mean(PLI_HC_theta,2),3),4),5);
mean_PLI_NO = mean(mean(mean(mean(PLI_NO_theta,2),3),4),5);

% Test average connectivity HC vs pt
[p,h,stats]=ranksum(mean_PLI_HC, mean_PLI_NO); % p=0.9075, U=3336 --> new: p=0.9356, U=3296

[Q,IQR] = quartile(mean_PLI_HC);
med = median(mean_PLI_HC);
[Q,IQR] = quartile(mean_PLI_NO);
med = median(mean_PLI_NO);

% Test relative power HC vs Pt
[p,h,stats] = ranksum(mean(mean_rel_pow_HC_theta_all,2), mean(mean_rel_pow_pt_theta,2)); % p=0.0537 U=2930 --> new: p=0.0690 U=2952

[Q,IQR] = quartile(mean(mean_rel_pow_HC_theta,2));
med = median(mean(mean_rel_pow_HC_theta,2));
[Q,IQR] = quartile(mean(mean_rel_pow_pt_theta,2));
med = median(mean(mean_rel_pow_pt_theta,2));

%% 4. Correlation theta power & global clus

mean_norm_clus_pt = squeeze(mean(mean(norm_cl_pt,3),2));

[rho, pval] = corr(mean_PLI_NO,mean_norm_clus_pt, 'type', 'Spearman'); % rho=0.0814, p=0.4992
[rho, pval] = corr(mean(mean_rel_pow_pt_theta,2), mean_norm_clus_pt, 'type', 'Spearman'); % rho = -0.0425 p=0.7246 

figure(); 
scatter(mean_PLI_NO, mean_norm_clus_pt)
lsline
%% 5. Tumor vs non-tumor pt

subj_pt = 1:71;
epochs_pt = 1:15; 
sub_epochs = 1:4;

for sub = 1:numel(subj_pt)
    for ep = 1:numel(epochs_pt)   
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(PLI_NO_theta(sub,ep,sub_ep,:,:));
            cl_pt_theta(sub,ep,sub_ep,:) = clustering_coef_wu(matrix);              % size matrix = 71x15x4x78 (clustering per region, per epoch, per patient)
        end
    end
end

% Norm-scores PLI & clustering patients
    
for sub = 1:numel(subj_pt)
    for ep = 1:numel(epochs_pt)
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(PLI_NO_theta(sub,ep,sub_ep,:,:));
            matrix_norm = weight_conversion(matrix, 'normalize');
            PLI_theta_norm_pt(sub,ep,sub_ep,:,:) = matrix_norm;                             % Normalized PLI matrix of each patient, size=71x15x4x78x78
            matrix_norm_PLI = squeeze(PLI_theta_norm_pt(sub,ep,sub_ep,:,:));
            clus_theta_norm_pt(sub,ep,sub_ep,:) = clustering_coef_wu(matrix_norm_PLI);        % Normalized local clus matrix of each patient, size=71x15x4x78
        end
    end
end

av_roi_cl_norm_pt = squeeze(mean(mean(clus_theta_norm_pt,3),2));                              % Mean normalized local clus per patient, averaged over epochs, size 71x78 

% Correlation average connectivity and clustering 

[rho, pval] = corr(mean_PLI_NO,mean(av_roi_cl_norm_pt,2), 'type', 'Spearman'); % rho=0.2985, p=0.0117
% figure(); 
% scatter(mean_PLI_NO, mean(av_roi_cl_norm_pt))
% lsline

% tumor roi and  non-tumor roi clustering (weighted for each tumor roi)

addpath('path/to/masks')                             % import data nr of voxels per AAL-MNI region

cd('path/to/masks/')                                
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
[p,h,stats]=signrank(tumor_clus_pt, non_tumor_clus_pt); % p=0.0425 W=1632

[Q,IQR] = quartile(tumor_clus_pt);
med = median(tumor_clus_pt);
[Q,IQR] = quartile(non_tumor_clus_pt);
med = median(non_tumor_clus_pt);

figure()
boxplot([tumor_clus_pt non_tumor_clus_pt],[zeros(1,71) ones(1,71)])

%% Mean and median # of tumor regions 

for k = 1:numel(subj)
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL2(:,:,k) = importdata(fullFileName);
    tumor_reg = find(voxels_AAL2(1:78,1,k) >0);
    nr_tumor_reg(k) = numel(tumor_reg);
end

mean_nr_tumor_reg = mean(nr_tumor_reg);
sd_nr_tumor_reg = std(nr_tumor_reg); 

%% 6. Z-score tumor vs non-tumor

subj_HC = 1:53;
epochs_HC = 1:13; 
    
% Norm scores HC

for sub = 1:numel(subj_HC)
    for ep = 1:numel(epochs_HC)   
        for sub_ep = 1:numel(sub_epochs)
        matrix = squeeze(PLI_HC_theta(sub,ep,sub_ep,:,:));
        matrix_norm = weight_conversion(matrix, 'normalize');
        PLI_theta_norm_HC(sub,ep,sub_ep,:,:) = matrix_norm;                     % Normalized PLI matrix of each HC, size=71x15x4x78x78
        matrix2 = squeeze(PLI_theta_norm_HC(sub,ep,sub_ep,:,:));
        clus_theta_norm_HC(sub,ep,sub_ep,:) = clustering_coef_wu(matrix2);      % Normalized local clus matrix of each HC, size=65x15x4x78
        end
    end
end

av_cl_norm_HC = squeeze(mean(mean(mean(clus_theta_norm_HC,4),3),2));            % Mean normalized local clus per HC, size=65x1
av_roi_cl_norm_HC = squeeze(mean(mean(mean(clus_theta_norm_HC))));              % Mean normalized local clus per roi, size=78x1
av_roi_cl_norm_HC2 = squeeze(mean(mean(mean(clus_theta_norm_HC,3),2),1));

% HC tumor roi and  non-tumor roi clustering (weighted for each tumor roi)

cd('path/to/masks/')

tumor_clus_HC = zeros(numel(subj),1);

for k = 1:numel(subj)  
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL = importdata(fullFileName);
    tumor_clus_HC(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*av_roi_cl_norm_HC);    % Clustering of weighted tumor regions taken together per HC, size 71x1
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    non_tumor_clus_HC(k,:) = mean(av_roi_cl_norm_HC(non_tumor_voxels_AAL));                     % Clustering of average of all non tumor regions per HC, size 71x1
end

av_cl_roi_HC=squeeze(mean(mean(clus_theta_norm_HC,3),2));                     % Local clustering for each region of each HC - 65x78
sd_av_cl_roi_HC = std(av_cl_roi_HC,0,1)';                                   % 78x1

for i = 1:71
    z_score(:,i) = (av_roi_cl_norm_pt(i,:)' - av_roi_cl_norm_HC)./sd_av_cl_roi_HC;  % clustering of each patient - clustering of HC regions/sd of hc regions
end 

z_score = z_score';

% z-score analyses

for k = 1:numel(subj) 
    fullFileName = fullfile(pwd, files(k).name);
    voxels_AAL = importdata(fullFileName);
    z_cl_roi_tumor(k,:) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*z_score(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    z_cl_roi_non_tumor(k,:) = mean(z_score(k,non_tumor_voxels_AAL));
end

% Test difference between tumor and non-tumor regions with z-scores
[p h stats] = signrank(z_cl_roi_tumor,z_cl_roi_non_tumor); % p=0.1537, W=1527 --> new: p=0.1138, W=1554

[Q,IQR] = quartile(z_cl_roi_tumor);
med = median(z_cl_roi_tumor);
[Q,IQR] = quartile(z_cl_roi_non_tumor);
med = median(z_cl_roi_non_tumor);

figure()
boxplot([z_cl_roi_tumor,z_cl_roi_non_tumor],[zeros(1,71) ones(1,71)]);
xlabel(['Intrinsic tumor clustering', 'Intrinsic non-tumor clustering'])
ylabel('Local clustering')

%% 7. Euclidian distance whole group

load('path/to/PLI_euclidian.mat')

res_euc_dis = reshape(euc_distance, 5538, [ ]);
res_hyp_dis = reshape(hyperbolic_distance,5538, [ ]);

z_score_flip = z_score';
res_cl_pt_roi_z = reshape(z_score_flip,5538,[ ]);

res_euc_dis_z(res_euc_dis<=0) = NaN;
res_hyp_dis_z(res_hyp_dis<=0) = NaN;

res_euc_dis(res_euc_dis<=0) = NaN;
res_hyp_dis(res_hyp_dis<=0) = NaN;

[rho_euc_z, pval_euc_z] = corr(res_euc_dis, res_cl_pt_roi_z,'rows','complete','type', 'Spearman'); % rho=0.0181, p=0.1770 --> new: rho=0.0184, p=0.1701
[rho_hyp_z, pval_hyp_z] = corr(res_hyp_dis, res_cl_pt_roi_z,'rows','complete','type', 'Spearman'); % rho=-.0127, p=0.3463 --> new: rho=-0.0120, p=0.3720

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
LME = fitlme(table(res_cl_pt_roi_z,res_euc_dis,grouping),'res_euc_dis~res_cl_pt_roi_z+(1|grouping)+(res_cl_pt_roi_z-1|grouping)'); %new p=0.36082, DF =5536, tstat=0.90025

% New figures relation distance & clustering
% Whole group
figure()
scatter(res_euc_dis, res_cl_pt_roi_z)
lsline 
h = lsline;
set(h(1),'color', 'r')

figure()
scatter(res_euc_dis, res_cl_pt_roi_z, 'r')
gscatter(res_cl_pt_roi_z,res_euc_dis,grouping)

    
figure()
scatter(res_euc_dis, res_cl_pt_roi_z, 'k'); hold on
scatter(res_euc_dis(1:78), res_cl_pt_roi_z(1:78), 'r', 'filled')
scatter(res_euc_dis(2*78-1:3*78), res_cl_pt_roi_z(2*78-1:3*78), 'g', 'filled')
scatter(res_euc_dis(5*78-1:6*78), res_cl_pt_roi_z(5*78-1:6*78), 'b', 'filled')

%% Figure LME

figure()
for i = 1:70
    scatter(res_euc_dis, res_cl_pt_roi_z, 'w'); hold on
    lsline; hold on
    scatter(res_euc_dis(i*78:(i+1)*78), res_cl_pt_roi_z(i*78:(i+1)*78), 'w')
    lsline
    hold on
%     h = lsline;
%     set(h(1),'color')
end

%% 8. Correlation tumor occurrence & local clus HC --> fig 2

load('path/to/PLI_overlap_N71.mat')

[rho, pval] = corr(av_roi_cl_norm_HC, aal_roi_tumor_overlap_N71, 'type', 'Spearman'); % p<0.001 en rho: 0.3821 --> new: p=0.0018, rho=0.3493

figure(); 
scatter(av_roi_cl_norm_HC, aal_roi_tumor_overlap_N71)
title('Tumor occurence related to local clustering of HC network')
xlabel('Average local clustering HC')
ylabel('Tumor occurence')
ylim([0 10])

%% 9. High vs low tumor z-scores    --> fig 3

mean_norm_glob_clus_pt = squeeze(mean(mean(norm_cl_pt,3),2));                % mean global clus per patient, size=71x1
norm_clus_pt_high = find(mean_norm_glob_clus_pt > 1.006881);                 % n=18
norm_clus_pt_low = find(mean_norm_glob_clus_pt < 1.006881);                  % n=53

[p h stats] = ranksum(z_cl_roi_tumor(norm_clus_pt_low),z_cl_roi_tumor(norm_clus_pt_high)); % p=0.0013, W=1664 --> new: p=0.0015, W=1667

[Q,IQR] = quartile(z_cl_roi_tumor(norm_clus_pt_low));
med = median(z_cl_roi_tumor(norm_clus_pt_low));
[Q,IQR] = quartile(z_cl_roi_tumor(norm_clus_pt_high));
med = median(z_cl_roi_tumor(norm_clus_pt_high));

figure()
boxplot([z_cl_roi_tumor(norm_clus_pt_low)' z_cl_roi_tumor(norm_clus_pt_high)'],[zeros(1,53) ones(1,18)])
xlabel(['tumor clustering low group', 'Tumor clustering high group'])
ylabel('local clustering value')

%% 10. High vs low non-tumor z-scores   --> fig 3 

[p h stats] = ranksum(z_cl_roi_non_tumor(norm_clus_pt_low), z_cl_roi_non_tumor(norm_clus_pt_high)); % p<0.001, W=1599

[Q,IQR] = quartile(z_cl_roi_non_tumor(norm_clus_pt_low));
med = median(z_cl_roi_non_tumor(norm_clus_pt_low));
[Q,IQR] = quartile(z_cl_roi_non_tumor(norm_clus_pt_high));
med = median(z_cl_roi_non_tumor(norm_clus_pt_high));

figure()
boxplot([z_cl_roi_non_tumor(norm_clus_pt_low)' z_cl_roi_non_tumor(norm_clus_pt_high)'],[zeros(1,53) ones(1,18)])
xlabel(['non-tumor clustering low group', 'non-tumor clustering high group'])
ylabel('local clustering value')

%% 11&12. High vs low relative power & average connectivity - tumor & non-tumor

cd('path/to/masks/')
files = dir('*_AAL.txt');

PLI_theta_mean = squeeze(mean(mean(mean(PLI_theta,5),3),2));

% Pt average connectivity n=18 
for k = 1:numel(norm_clus_pt_high)  
    fullFileName = fullfile(pwd, files(norm_clus_pt_high(k)).name);
    voxels_AAL = importdata(fullFileName);
    PLI_theta_pt_tumor18(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*PLI_theta_mean(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    PLI_theta_pt_non_tumor18(k,:) = mean(PLI_theta_mean(k,non_tumor_voxels_AAL));
end

% Pt average connectivity n=53 
for k = 1:numel(norm_clus_pt_low)  
    fullFileName = fullfile(pwd, files(norm_clus_pt_low(k)).name);
    voxels_AAL = importdata(fullFileName);
    PLI_theta_pt_tumor53(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*PLI_theta_mean(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    PLI_theta_pt_non_tumor53(k,:) = mean(PLI_theta_mean(k,non_tumor_voxels_AAL));
end

% Pt relative power n=18
for k = 1:numel(norm_clus_pt_high)  
    fullFileName = fullfile(pwd, files(norm_clus_pt_high(k)).name);
    voxels_AAL = importdata(fullFileName);
    rel_pow_theta_pt_tumor18(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*mean_rel_pow_pt_theta(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    rel_pow_theta_pt_non_tumor18(k,:) = mean(mean_rel_pow_pt_theta(k,non_tumor_voxels_AAL));
end

% Pt relative power n=53
for k = 1:numel(norm_clus_pt_low)  
    fullFileName = fullfile(pwd, files(norm_clus_pt_low(k)).name);
    voxels_AAL = importdata(fullFileName);
    rel_pow_theta_pt_tumor53(k,1) = sum(voxels_AAL(1:78,1)/sum(voxels_AAL(1:78,1)).*mean_rel_pow_pt_theta(k,:)');
    non_tumor_voxels_AAL = find(voxels_AAL(1:78,1)==0);
    rel_pow_theta_pt_non_tumor53(k,:) = mean(mean_rel_pow_pt_theta(k,non_tumor_voxels_AAL));
end

% Test tumor & non-tumor n=18 vs n=53 - power 
[p,h,stats]=ranksum(rel_pow_theta_pt_tumor18, rel_pow_theta_pt_tumor53); % p=0.9210 U=656 --> zelfde
[p,h,stats]=ranksum(rel_pow_theta_pt_non_tumor18, rel_pow_theta_pt_non_tumor53); % p=0.4632 U=704 --> zelfde

[Q,IQR] = quartile(rel_pow_theta_pt_tumor18);
med = median(rel_pow_theta_pt_tumor18);
[Q,IQR] = quartile(rel_pow_theta_pt_tumor53);
med = median(rel_pow_theta_pt_tumor53);

[Q,IQR] = quartile(rel_pow_theta_pt_non_tumor18);
med = median(rel_pow_theta_pt_non_tumor18);
[Q,IQR] = quartile(rel_pow_theta_pt_non_tumor53);
med = median(rel_pow_theta_pt_non_tumor53);

% Test tumor & non-tumor n=18 vs n=53 - average connectivity
[p,h,stats]=ranksum(PLI_theta_pt_tumor18, PLI_theta_pt_tumor53); % p = 0.1592 U=755 
[p,h,stats]=ranksum(PLI_theta_pt_non_tumor18, PLI_theta_pt_non_tumor53); % p = 0.3866 U = 582 

[Q,IQR] = quartile(PLI_theta_pt_tumor18);
med = median(PLI_theta_pt_tumor18);
[Q,IQR] = quartile(PLI_theta_pt_tumor53);
med = median(PLI_theta_pt_tumor53);

[Q,IQR] = quartile(PLI_theta_pt_non_tumor18);
med = median(PLI_theta_pt_non_tumor18);
[Q,IQR] = quartile(PLI_theta_pt_non_tumor53);
med = median(PLI_theta_pt_non_tumor53);


%% 13. High vs low HC tumor --> fig 4

% Test difference between tumor and non tumor clustering in HC network
[p,h,stats]=ranksum(tumor_clus_HC(norm_clus_pt_high), tumor_clus_HC(norm_clus_pt_low)); % p=0.0170 U=467 --> new: p=0.0579, U=504

[Q,IQR] = quartile(tumor_clus_HC(norm_clus_pt_high));
med = median(tumor_clus_HC(norm_clus_pt_high));
[Q,IQR] = quartile(tumor_clus_HC(norm_clus_pt_low));
med = median(tumor_clus_HC(norm_clus_pt_low));

figure()
boxplot([tumor_clus_HC(norm_clus_pt_high)', tumor_clus_HC(norm_clus_pt_low)'],[zeros(1,18) ones(1,53)])
xlabel(['tumor clustering HC high clustering group', 'tumor clustering HC low clustering group'])
ylabel('local clustering value')

%% 14. Euclidian distance high & low

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

[rho_euc_z_high, pval_euc_z_high] = corr(res_euc_dis_tumor_high, res_z_score_high, 'type', 'Spearman'); % rho= 0.0763, p=0.0042 --> new: p=0.0044, rho=0.0759
[rho_euc_z_low, pval_euc_z_low] = corr(res_euc_dis_tumor_low, res_z_score_low, 'type', 'Spearman'); % rho=0.006 , p=0.6988 --> new: p=0.7662. rho=0.0046
[rho_hyp_z_high, pval_hyp_z_high] = corr(res_hyp_dis_tumor_high, res_z_score_high, 'type', 'Spearman'); % rho= -0.0389 , p= 0.1447 --> p=0.1262, rho=-0.0408
[rho_hyp_z_low, pval_hyp_z_low] = corr(res_hyp_dis_tumor_low, res_z_score_low, 'type', 'Spearman'); % rho= -0.0197, p=0.2057 --> new: p=0.1871, rho=-0.0205

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

% New figures relation distance & clustering
% high group
figure()
scatter(res_euc_dis_tumor_high,res_z_score_high)
lsline 
h = lsline;
set(h(1),'color', 'r')
ylim([-4 6])

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

sig_euc_dis_ind=find(pval_euc_z_ind <0.08);

for i=1:20
    figure(); 
    scatter(res_euc_dis_ind(:,sig_euc_dis_ind(i)),res_cl_pt_roi_z_ind(:,sig_euc_dis_ind(i)))
    title('Individual euclidian distance in relation to local clustering')
    lsline
end

%% Relationship volume with global and tumor clustering

tumor_volume = [48174 40055 57008 16559 5284 13915 86604 4764 123054 869 57557 46525 913 111484 66194 9118 9641 11545 5961 76151 14226 3378 62942 12785 73978 8863 5617 8213 13493 36209 6815 11088 18569 17016 17922 3538 43732 62804 53097 58605 69134 51568 40576 55132 31840 38992 48471 47304 17651 31217 22187 25495 32859 127000 42615 11013 36761 41015 67971 62346 99097 102514 1600 107067 15062 30097 19402 23027 66208 212018 25284];
%From SPSS, same order as matlab matrix 

[rho_cor_glob_vol p_cor_glob_vol] = corr(norm_cl_NO,tumor_volume', 'Type', 'Spearman');
[rho_cor_tum_vol p_cor_tum_vol] = corr(tumor_clus_pt,tumor_volume', 'Type', 'Spearman');

figure; 
scatter(norm_cl_NO,tumor_volume')
xlabel('Global clustering glioma patients')
ylabel('Tumor volume')

figure; 
scatter(tumor_clus_pt,tumor_volume')
xlabel('Local tumor clustering')
ylabel('Tumor volume')

load('path/to/PLI_rel_tum_vol.mat')
delete_vs = [2 43 56 66]; 
norm_cl_NO(delete_vs) = []; 
tumor_clus_pt(delete_vs) = []; 

[rho_cor_glob_rel_vol p_cor_glob_rel_vol] = corr(norm_cl_NO,rel_tum_vol', 'Type', 'Spearman');
[rho_cor_tum_rel_vol p_cor_tum_rel_vol] = corr(tumor_clus_pt,rel_tum_vol', 'Type', 'Spearman');

figure; 
scatter(norm_cl_NO,rel_tum_vol')
lsline
xlabel('Global clustering glioma patients')
ylabel('Relative tumor volume')

figure; 
scatter(tumor_clus_pt,rel_tum_vol')
lsline
xlabel('Local tumor clustering')
ylabel('Relative tumor volume')
