% Multiscale integration and memory in TLE

%   This script is used to do network analysis on both fMRI and MEG time
%   series from a cohort of TLE patients, and combines these analyses with
%   available micro-scale and memory variables.
%   The paper accompanying this analysis can be found at
%   https://academic.oup.com/cercor/advance-article/doi/10.1093/cercor/bhab349/6375261 

%%%

%%author__ = Linda Douw
%%contact__ = l.douw@amsterdamumc.nl
%%date__ = 2021/10/04
%%status__ = finished 

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

%%% Reviewed by Name Date 

% NA

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% NA

%%% Toolboxes 
% Brain Connectivity Toolbox https://sites.google.com/site/bctnet/ 
% version 2017-15-01

%% add local path for BCT and kruskal algorithm

% addpath /path/bct
% addpath /path/x % contains kruskal_algorithm.m
% 
% %% Read in adjacency matrices for fMRI and MEG, clinical and micro data, 
% % and data availability overview
% 
% load /path/aal_fmri_adjmats.mat
% load /path/aal_meg_adjmats_t.mat
% load /path/aal_func_multilayer_ecm.mat
% load /path/micro_clinical_data.mat
% load /path/additional_analyses.mat

% Note that non-anonymizable entries in the clinical/micro datafiles are
% censored using NA (for variable names) and NaN (for elements)

%% Create fMRI-related variables

% Number of regions and subjects
rois = 72;
nsubs_fmri = 20; 

% Var indicating lateralization of the TLE and resection
lat_fmri = clin_micro_data(find(data_availability_overview(:,2)),9);

% Extract network signatures of the resected area (middle temporal gyrus) 
% AAL region 64 in the right hemisphere, region 28 in the left hemisphere
aal_regions_fmri = lat_fmri;
aal_regions_fmri(lat_fmri == 1,:) = 28;
aal_regions_fmri(lat_fmri == 2,:) = 64;

%% Calculate MST on fMRI adjacency matrices

aal_fmri_adjmats_mst = zeros(rois,rois,nsubs_fmri);
for i = 1:nsubs_fmri
    temp = 1./aal_fmri_adjmats(:,:,i); % minimum instead of maximum spanning tree
    temp(isinf(temp)) = 0;
    [tempmat,~,~] = kruskal_algorithm(temp);
    aal_fmri_adjmats_mst(:,:,i) = tempmat;
end

% calculate EC per node and specifically for the resected area
aal_fmri_ecm_mst = nan(nsubs_fmri,rois);
aal_fmri_ecm_mst_resar = nan(nsubs_fmri,1);
for i = 1:nsubs_fmri 
    aal_fmri_ecm_mst(i,:) = (eigenvector_centrality_und(aal_fmri_adjmats_mst(:,:,i)))';
    aal_fmri_ecm_mst_resar(i,:) = aal_fmri_ecm_mst(i,aal_regions_fmri(i,:));
end

%% Select subsamples for fMRI x TDL / AP rise speeds correlations

subs_fmri_tdl = nan(length(data_availability_overview),3);
subs_fmri_tdl(:,1) = data_availability_overview(:,2);
subs_fmri_tdl(:,2) = data_availability_overview(:,6);
subs_fmri_tdl(:,3) = sum(data_availability_overview(:,[2 6]),2);

% combine fMRI and TDL 
select_subs_tdl_for_fmri = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_fmri_tdl(i,2) == 1 && subs_fmri_tdl(i,3) == 2
        select_subs_tdl_for_fmri(i,:) = 1;
    end
end

temp_select_subs_fmri_for_tdl = find(data_availability_overview(:,2));
select_subs_fmri_for_tdl = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_tdl(temp_select_subs_fmri_for_tdl(i,:),2) == 1
        select_subs_fmri_for_tdl(i,:) = 1;
    end
end

clear temp* i

subs_fmri_ephys = nan(length(data_availability_overview),3);
subs_fmri_ephys(:,1) = data_availability_overview(:,2);
subs_fmri_ephys(:,2) = data_availability_overview(:,5);
subs_fmri_ephys(:,3) = sum(data_availability_overview(:,[2 5]),2);

select_subs_ephys_for_fmri = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_fmri_ephys(i,2) == 1 && subs_fmri_ephys(i,3) == 2
        select_subs_ephys_for_fmri(i,:) = 1;
    end
end

temp_select_subs_fmri_for_ephys = find(data_availability_overview(:,2));
select_subs_fmri_for_ephys = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_ephys(temp_select_subs_fmri_for_ephys(i,:),2) == 1
        select_subs_fmri_for_ephys(i,:) = 1;
    end
end

clear temp* i

%% Perform correlation of fMRI EC with TDL / AP rise speeds

[corr_fmri_ecm_mst_tdl_resar,corr_fmri_ecm_mst_tdl_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_fmri),10),'type','spearman');
[corr_fmri_ecm_mst_aprise1st_resar,corr_fmri_ecm_mst_aprise1st_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),2),'type','spearman');
[corr_fmri_ecm_mst_aprise2040_resar,corr_fmri_ecm_mst_aprise2040_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),12),'type','spearman');

%% Validation of significant results
% First explore regional specificity of significant correlations to resected node

aal_fmri_adjmats_mst_ipsi = nan(rois/2,rois,nsubs_fmri);
aal_fmri_ecm_mst_ipsi = nan(nsubs_fmri,rois/2);
for i = 1:nsubs_fmri
    if lat_fmri(i,:) == 1
        aal_fmri_ecm_mst_ipsi(i,:) = aal_fmri_ecm_mst(i,1:rois/2);
        aal_fmri_adjmats_mst_ipsi(:,:,i) = aal_fmri_adjmats_mst(1:rois/2,:,i);
    end
    if lat_fmri(i,:) == 2
        aal_fmri_ecm_mst_ipsi(i,:) = aal_fmri_ecm_mst(i,rois/2+1:end);
        aal_fmri_adjmats_mst_ipsi(:,:,i) = aal_fmri_adjmats_mst(rois/2+1:end,:,i);
    end
end

corr_fmri_ecm_mst_tdl_ipsi = nan(rois/2,2);
corr_fmri_ecm_mst_aprise1st_ipsi = nan(rois/2,2);

for i = 1:rois/2    
    [corr_fmri_ecm_mst_tdl_ipsi(i,1),corr_fmri_ecm_mst_tdl_ipsi(i,2)] = ...
        corr(aal_fmri_ecm_mst_ipsi(find(select_subs_fmri_for_tdl),i),...
        clin_micro_data(find(select_subs_tdl_for_fmri),10), ...
        'type','spearman');
    [corr_fmri_ecm_mst_aprise1st_ipsi(i,1),corr_fmri_ecm_mst_aprise1st_ipsi(i,2)] = ...
        corr(aal_fmri_ecm_mst_ipsi(find(select_subs_fmri_for_ephys),i),...
        clin_micro_data(find(select_subs_ephys_for_fmri),2), ...
        'type','spearman');
end

% Create histogram of significant correlations (figure 3B / 3D)

hist(corr_fmri_ecm_mst_tdl_ipsi(:,1),15,'FaceColor','none')
box off

hist(corr_fmri_ecm_mst_aprise1st_ipsi(:,1),15,'FaceColor','none')
box off

%% Check significant correlations against sample-specific distribution through permutation (figure 3B/3D)

% for TDL 
temp_aal_fmri_ecm_mst = aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_tdl),:);
temp_tdl = clin_micro_data(find(select_subs_tdl_for_fmri),10);

nperms = 10000;

corrs_perm_aal_fmri_ecm_mst_tdl_resar = nan(nperms,2);
for i = 1:nperms
    temp_perm_order = randperm(10);
    [corrs_perm_aal_fmri_ecm_mst_tdl_resar(i,1),corrs_perm_aal_fmri_ecm_mst_tdl_resar(i,2)] = ...
        corr(temp_aal_fmri_ecm_mst(temp_perm_order),temp_tdl,'type','spearman');
end

temp_perc = prctile(corrs_perm_aal_fmri_ecm_mst_tdl_resar(:,1),97.5);

hist(corrs_perm_aal_fmri_ecm_mst_tdl_resar(:,1),50)
hold on;
line([0.758, 0.758], ylim, 'LineWidth', 2, 'Color', 'g'); % real correlation
hold on;
line([0.6364, 0.6364], ylim, 'LineWidth', 2, 'Color', 'r'); % temp_perc cut-off
xlabel('permutated correlations (n=10000) between fMRI ECM and TDL in mm')
ylabel('frequency of occurrence')

clear temp* i

% For AP rise speed

temp_aal_fmri_ecm = aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:);
temp_ephys = clin_micro_data(find(select_subs_ephys_for_fmri),2);

nperms = 10000;

corrs_perm_aal_fmri_ecm_mst_aprise1st_resar = nan(nperms,2);
for i = 1:nperms
    temp_perm_order = randperm(15);
    [corrs_perm_aal_fmri_ecm_mst_aprise1st_resar(i,1),corrs_perm_aal_fmri_ecm_mst_aprise1st_resar(i,2)] = ...
        corr(temp_aal_fmri_ecm(temp_perm_order),temp_ephys,'type','spearman');
end

temp_perc = prctile(corrs_perm_aal_fmri_ecm_mst_aprise1st_resar(:,1),97.5);

hist(corrs_perm_aal_fmri_ecm_mst_aprise1st_resar(:,1),50)
hold on;
line([0.539, 0.539], ylim, 'LineWidth', 2, 'Color', 'g');
hold on;
line([0.5268, 0.5268], ylim, 'LineWidth', 2, 'Color', 'r');

clear temp* i

%% Leave-one-out-validation of significant results (figure 3B/3D)

corrs_loo_aal_fmri_ecm_mst_tdl_resar = nan(sum(select_subs_tdl_for_fmri),2);
for i = 1:sum(select_subs_tdl_for_fmri)
    temp_aal_fmri_ecm = aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_tdl),:);
    temp_tdl = clin_micro_data(find(select_subs_tdl_for_fmri),10);
    temp_aal_fmri_ecm(i) = [];
    temp_tdl(i) = [];
    [corrs_loo_aal_fmri_ecm_mst_tdl_resar(i,1),corrs_loo_aal_fmri_ecm_mst_tdl_resar(i,2)] = ...
        corr(temp_aal_fmri_ecm,temp_tdl,'type','spearman');
end

corrs_loo_aal_fmri_ecm_mst_aprise1st_resar = nan(sum(select_subs_ephys_for_fmri),2);
for i = 1:sum(select_subs_ephys_for_fmri)
    temp_aal_fmri_ecm = aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:);
    temp_aprise1st = clin_micro_data(find(select_subs_ephys_for_fmri),2);
    temp_aal_fmri_ecm(i) = [];
    temp_aprise1st(i) = [];
    [corrs_loo_aal_fmri_ecm_mst_aprise1st_resar(i,1),corrs_loo_aal_fmri_ecm_mst_aprise1st_resar(i,2)] = ...
        corr(temp_aal_fmri_ecm,temp_aprise1st,'type','spearman');
end

clear temp* i 

hist(corrs_loo_aal_fmri_ecm_mst_tdl_resar(:,1),10)
yticks([0 1 2 3 4])
box off 

figure
hist(corrs_loo_aal_fmri_ecm_mst_aprise1st_resar(:,1),10)
%yticks([0 1 2 3 4])
box off 

%% MEG analysis
% First create MEG variables 
% Since there are 78 (not 72, like with fMRI) regions in the MEG data:
% AAL region 70 in the right hemisphere, region 31 in the left hemisphere

nsubs_meg = sum(data_availability_overview(:,3));
rois = 78;

lat_meg = clin_micro_data(find(data_availability_overview(:,3)),9);
aal_regions_meg = nan(nsubs_meg,1);
for i = 1:nsubs_meg
    if lat_meg(i,:) == 1 
        aal_regions_meg(i,:) = 31;
    elseif lat_meg(i,:) == 2 
        aal_regions_meg(i,:) = 70;
    end
end

%% Select subsamples for MEG x TDL / AP rise speed correlations

subs_meg_tdl = nan(length(data_availability_overview),3);
subs_meg_tdl(:,1) = data_availability_overview(:,3);
subs_meg_tdl(:,2) = data_availability_overview(:,6);
subs_meg_tdl(:,3) = sum(data_availability_overview(:,[3 6]),2);

select_subs_tdl_for_meg = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_meg_tdl(i,2) == 1 && subs_meg_tdl(i,3) == 2
        select_subs_tdl_for_meg(i,:) = 1;
    end
end

temp_select_subs_meg_for_tdl = find(data_availability_overview(:,3));
select_subs_meg_for_tdl = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_tdl(temp_select_subs_meg_for_tdl(i,:),2) == 1
        select_subs_meg_for_tdl(i,:) = 1;
    end
end

clear temp* i

subs_meg_ephys = nan(length(data_availability_overview),3);
subs_meg_ephys(:,1) = data_availability_overview(:,3);
subs_meg_ephys(:,2) = data_availability_overview(:,5);
subs_meg_ephys(:,3) = sum(data_availability_overview(:,[3 5]),2);

select_subs_ephys_for_meg = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_meg_ephys(i,2) == 1 && subs_meg_ephys(i,3) == 2
        select_subs_ephys_for_meg(i,:) = 1;
    end
end

temp_select_subs_meg_for_ephys = find(data_availability_overview(:,3));
select_subs_meg_for_ephys = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_ephys(temp_select_subs_meg_for_ephys(i,:),2) == 1
        select_subs_meg_for_ephys(i,:) = 1;
    end
end

clear temp* i

%% Calculate MST on MEG adjacency matrices

aal_meg_adjmats_t_mst = zeros(rois,rois,nsubs_meg);
for i = 1:nsubs_meg
    temp = 1./aal_meg_adjmats_t(:,:,i); % minimum instead of maximum spanning tree
    temp(isinf(temp)) = 0;
    [tempmat,~,~] = kruskal_algorithm(temp);
    aal_meg_adjmats_t_mst(:,:,i) = tempmat;
end

%% Calculate MST EC per node and specifically for resected region

aal_meg_t_ecm_mst = nan(nsubs_meg,rois);
aal_meg_t_ecm_mst_resar = nan(nsubs_meg,1);
for i = 1:nsubs_meg 
    aal_meg_t_ecm_mst(i,:) = eigenvector_centrality_und(aal_meg_adjmats_t_mst(:,:,i));
    aal_meg_t_ecm_mst_resar(i,:) = aal_meg_t_ecm_mst(i,aal_regions_meg(i,:));
end

%% Correlate MEG EC with micro-scale TDL / AP rise speeds

[corr_meg_t_ecm_mst_tdl_resar,corr_meg_t_ecm_mst_tdl_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');
[corr_meg_t_ecm_mst_aprise1st_resar,corr_meg_t_ecm_mst_aprise1st_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_meg_t_ecm_mst_aprise2040_resar,corr_meg_t_ecm_mst_aprise2040_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');

%% We here create the multilayer supra adjacency matrix
% First select sample with both fMRI and MEG

rois = 72; % because only 72 regions in fMRI

subs_func_multilayer = nan(length(data_availability_overview),3);
subs_func_multilayer(:,1) = data_availability_overview(:,2);
subs_func_multilayer(:,2) = data_availability_overview(:,3);
subs_func_multilayer(:,3) = sum(data_availability_overview(:,2:3),2);

temp_select_subs_fmri_for_func_multilayer = find(data_availability_overview(:,2));
select_subs_fmri_for_func_multilayer = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_func_multilayer(temp_select_subs_fmri_for_func_multilayer(i,:),3) == 2
        select_subs_fmri_for_func_multilayer(i,:) = 1;
    end
end

nsubs_func_multilayer = sum(select_subs_fmri_for_func_multilayer);

temp_select_subs_meg_for_func_multilayer = find(data_availability_overview(:,3));
select_subs_meg_for_func_multilayer = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_func_multilayer(temp_select_subs_meg_for_func_multilayer(i,:),3) == 2
        select_subs_meg_for_func_multilayer(i,:) = 1;
    end
end

clear temp* i

%% Save these supra adjacency matrices for analysis in python

%first remove missing fMRI regions from meg
aal_fmri_regions_nosignal = [1,2,4,40,41,42]; % this is the result from extensive visual inspection
aal_meg_adjmats_t_mst_72 = aal_meg_adjmats_t_mst;
aal_meg_adjmats_t_mst_72(aal_fmri_regions_nosignal,:,:) = [];
aal_meg_adjmats_t_mst_72(:,aal_fmri_regions_nosignal,:) = [];

tempdiag = ones(rois,1);
tempmat = diag(tempdiag);

temp_fmri_for_func_multilayer = aal_fmri_adjmats_mst(:,:,select_subs_fmri_for_func_multilayer == 1);
temp_meg_for_func_multilayer = aal_meg_adjmats_t_mst_72(:,:,select_subs_meg_for_func_multilayer == 1);

aal_supra_adjmat_mst = nan(rois*2,rois*2,nsubs_func_multilayer);
for i = 1:nsubs_func_multilayer
    aal_supra_adjmat_mst(1:rois,1:rois,i) = temp_fmri_for_func_multilayer(:,:,i);
    aal_supra_adjmat_mst(rois+1:end,rois+1:end,i) = temp_meg_for_func_multilayer(:,:,i);
    aal_supra_adjmat_mst(1:rois,rois+1:end,i) = tempmat;
    aal_supra_adjmat_mst(rois+1:end,1:rois,i) = tempmat;
end

% aal_supra_adjmat_mst was used as input in our python multilayer pipeline,
% also on this GitHub

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO ANALYSIS OF MULTILAYER ECM IN PYTHON 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract multilayer EC of resected region
% remember: AAL region 64 in the right hemisphere, region 28 in the left hemisphere

lat_func_multilayer = clin_micro_data(find(data_availability_overview(:,11)),9);

aal_regions_func_multilayer = nan(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if lat_func_multilayer(i,:) == 1 
        aal_regions_func_multilayer(i,:) = 28;
    elseif lat_func_multilayer(i,:) == 2 
        aal_regions_func_multilayer(i,:) = 64;
    end
end

aal_func_multilayer_ecm_resar = nan(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    aal_func_multilayer_ecm_resar(i,:) = aal_func_multilayer_ecm(aal_regions_func_multilayer(i,:),i);
end

%% Select subsamples for multilayer x TDL / rise speeds correlations

subs_func_multilayer_tdl = nan(length(data_availability_overview),3);
subs_func_multilayer_tdl(:,1) = data_availability_overview(:,11); % multilayer
subs_func_multilayer_tdl(:,2) = data_availability_overview(:,6); % tdl
subs_func_multilayer_tdl(:,3) = sum(data_availability_overview(:,[11 6]),2);

select_subs_tdl_for_func_multilayer = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_func_multilayer_tdl(i,2) == 1 && subs_func_multilayer_tdl(i,3) == 2
        select_subs_tdl_for_func_multilayer(i,:) = 1;
    end
end

temp_select_subs_func_multilayer_for_tdl = find(data_availability_overview(:,11));
select_subs_func_multilayer_for_tdl = zeros(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if subs_func_multilayer_tdl(temp_select_subs_func_multilayer_for_tdl(i,:),2) == 1
        select_subs_func_multilayer_for_tdl(i,:) = 1;
    end
end

subs_func_multilayer_ephys = nan(length(data_availability_overview),3);
subs_func_multilayer_ephys(:,1) = data_availability_overview(:,11); % multilayer
subs_func_multilayer_ephys(:,2) = data_availability_overview(:,5); % ephys
subs_func_multilayer_ephys(:,3) = sum(data_availability_overview(:,[11 5]),2);

select_subs_ephys_for_func_multilayer = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_func_multilayer_ephys(i,2) == 1 && subs_func_multilayer_ephys(i,3) == 2
        select_subs_ephys_for_func_multilayer(i,:) = 1;
    end
end

temp_select_subs_func_multilayer_for_ephys = find(data_availability_overview(:,11));
select_subs_func_multilayer_for_ephys = zeros(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if subs_func_multilayer_ephys(temp_select_subs_func_multilayer_for_ephys(i,:),2) == 1
        select_subs_func_multilayer_for_ephys(i,:) = 1;
    end
end

clear temp*

%% Correlations 

[corr_func_multilayer_ecm_mst_aprise1st_resar,corr_func_multilayer_ecm_mst_aprise1st_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),2),...
    'type','spearman');
[corr_func_multilayer_ecm_mst_aprise2040_resar,corr_func_multilayer_ecm_mst_aprise2040_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),12),...
    'type','spearman');
[corr_func_multilayer_ecm_tdl_resar,corr_func_multilayer_ecm_tdl_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_tdl)), ...
    clin_micro_data(find(select_subs_tdl_for_func_multilayer),10),...
    'type','spearman');

%% Explore regional specificity of significant correlation to resected node (figure 4B)

aal_func_multilayer_ecm_ipsi = nan(nsubs_func_multilayer,rois/2);
for i = 1:nsubs_func_multilayer
    if lat_func_multilayer(i,:) == 1
        aal_func_multilayer_ecm_ipsi(i,:) = aal_func_multilayer_ecm(1:rois/2,i);
    elseif lat_func_multilayer(i,:) == 2
        aal_func_multilayer_ecm_ipsi(i,:) = aal_func_multilayer_ecm((rois/2)+1:end,i);
    end
end

corr_func_multilayer_ecm_mst_aprise1st_ipsi = nan(rois/2,2);
for i = 1:rois/2    
    [corr_func_multilayer_ecm_mst_aprise1st_ipsi(i,1),corr_func_multilayer_ecm_mst_aprise1st_ipsi(i,2)] = ...
        corr(aal_func_multilayer_ecm_ipsi(find(select_subs_func_multilayer_for_ephys),i),...
        clin_micro_data(find(select_subs_ephys_for_func_multilayer),2),'type','spearman');
end

hist(corr_func_multilayer_ecm_mst_aprise1st_ipsi(:,1),15)
yticks([0 1 2 3 4 5 6 7])
box off

%% Check significant correlations against sample-specific distribution through permutation (figure 4B)

temp_aal_func_multilayer_ecm = aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys),:);
temp_aprise1st = clin_micro_data(find(select_subs_ephys_for_func_multilayer),2);

nperms = 10000;

corrs_perm_aal_func_multilayer_ecm_aprise1st_resar = nan(nperms,2);
for i = 1:nperms
    temp_perm_order = randperm(7);
    [corrs_perm_aal_func_multilayer_ecm_aprise1st_resar(i,1),corrs_perm_aal_func_multilayer_ecm_aprise1st_resar(i,2)] = ...
        corr(temp_aal_func_multilayer_ecm(temp_perm_order),temp_aprise1st,'type','spearman');
end

temp_perc = prctile(corrs_perm_aal_func_multilayer_ecm_aprise1st_resar(:,1),97.5);

hist(corrs_perm_aal_func_multilayer_ecm_aprise1st_resar(:,1),50)
box off

clear temp* i
 
%% Leave-one-out-validation (figure 4B)

corrs_loo_aal_func_multilayer_ecm_resar = nan(sum(select_subs_func_multilayer_for_ephys),2);
for i = 1:sum(select_subs_func_multilayer_for_ephys)
    temp_aal_func_multilayer_ecm_resar = aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys),:);
    temp_aprise1st = clin_micro_data(find(select_subs_ephys_for_func_multilayer),2);
    temp_aal_func_multilayer_ecm_resar(i) = [];
    temp_aprise1st(i) = [];
    [corrs_loo_aal_func_multilayer_ecm_resar(i,1),corrs_loo_aal_func_multilayer_ecm_resar(i,2)] = ...
        corr(temp_aal_func_multilayer_ecm_resar,temp_aprise1st,'type','spearman');
end

hist(corrs_loo_aal_func_multilayer_ecm_resar(:,1),10)
box off
yticks([0 1 2 3 4 5])
xlim([0.90 1.03])
xticks([0.9 0.95 1])

%% Re-analyze unilayer correlations within multilayer sample only
% Select smaller samples from fMRI and MEG samples

subs_fmri_func_multilayer_ephys = nan(length(data_availability_overview),3);
subs_fmri_func_multilayer_ephys(:,1) = data_availability_overview(:,11);
subs_fmri_func_multilayer_ephys(:,2) = data_availability_overview(:,2); 
subs_fmri_func_multilayer_ephys(:,3) = data_availability_overview(:,5); 
subs_fmri_func_multilayer_ephys(:,4) = sum(data_availability_overview(:,[11 2 5]),2);

temp_select_subs_fmri_for_ephys_func_multilayer = find(data_availability_overview(:,2));
select_subs_fmri_for_ephys_func_multilayer = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_func_multilayer_ephys(temp_select_subs_fmri_for_ephys_func_multilayer(i,:),1) == 1 && ...
            subs_fmri_func_multilayer_ephys(temp_select_subs_fmri_for_ephys_func_multilayer(i,:),4) == 3
        select_subs_fmri_for_ephys_func_multilayer(i,:) = 1;
    end
end

subs_fmri_func_multilayer_tdl = nan(length(data_availability_overview),3);
subs_fmri_func_multilayer_tdl(:,1) = data_availability_overview(:,11);
subs_fmri_func_multilayer_tdl(:,2) = data_availability_overview(:,2); 
subs_fmri_func_multilayer_tdl(:,3) = data_availability_overview(:,6); 
subs_fmri_func_multilayer_tdl(:,4) = sum(data_availability_overview(:,[11 2 6]),2);

temp_select_subs_fmri_for_tdl_func_multilayer = find(data_availability_overview(:,2));
select_subs_fmri_for_tdl_func_multilayer = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_func_multilayer_tdl(temp_select_subs_fmri_for_tdl_func_multilayer(i,:),1) == 1 && ...
            subs_fmri_func_multilayer_tdl(temp_select_subs_fmri_for_tdl_func_multilayer(i,:),4) == 3
        select_subs_fmri_for_tdl_func_multilayer(i,:) = 1;
    end
end

subs_meg_func_multilayer_tdl = nan(length(data_availability_overview),3);
subs_meg_func_multilayer_tdl(:,1) = data_availability_overview(:,11); % multilayer
subs_meg_func_multilayer_tdl(:,2) = data_availability_overview(:,3); % meg
subs_meg_func_multilayer_tdl(:,3) = data_availability_overview(:,6); % tdl
subs_meg_func_multilayer_tdl(:,4) = sum(data_availability_overview(:,[11 3 6]),2);

temp_select_subs_meg_for_tdl_func_multilayer = find(data_availability_overview(:,3));
select_subs_meg_for_tdl_func_multilayer = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_func_multilayer_tdl(temp_select_subs_meg_for_tdl_func_multilayer(i,:),1) == 1 && ...
            subs_meg_func_multilayer_tdl(temp_select_subs_meg_for_tdl_func_multilayer(i,:),3) == 1
        select_subs_meg_for_tdl_func_multilayer(i,:) = 1;
    end
end

subs_meg_func_multilayer_ephys = nan(length(data_availability_overview),3);
subs_meg_func_multilayer_ephys(:,1) = data_availability_overview(:,11); % multilayer
subs_meg_func_multilayer_ephys(:,2) = data_availability_overview(:,3); % meg
subs_meg_func_multilayer_ephys(:,3) = data_availability_overview(:,5); % ephys
subs_meg_func_multilayer_ephys(:,4) = sum(data_availability_overview(:,[11 3 5]),2);

temp_select_subs_meg_for_ephys_func_multilayer = find(data_availability_overview(:,3));
select_subs_meg_for_ephys_func_multilayer = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_func_multilayer_ephys(temp_select_subs_meg_for_ephys_func_multilayer(i,:),1) == 1 && ...
            subs_meg_func_multilayer_ephys(temp_select_subs_meg_for_ephys_func_multilayer(i,:),4) == 3
        select_subs_meg_for_ephys_func_multilayer(i,:) = 1;
    end
end

%% Correlations in subsamples 

% fMRI
[corr_fmri_ecm_mst_tdl_resar_subgroup_fm,corr_fmri_ecm_mst_tdl_resar_subgroup_fm_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_tdl_func_multilayer)), ...
    clin_micro_data(find(select_subs_tdl_for_func_multilayer),10),...
    'type','spearman');
[corr_fmri_ecm_mst_aprise2040_resar_subgroup_fm,corr_fmri_ecm_mst_aprise2040_resar_subgroup_fm_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys_func_multilayer)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),12),...
    'type','spearman');
[corr_fmri_ecm_mst_aprise1st_resar_subgroup_fm,corr_fmri_ecm_mst_aprise1st_resar_subgroup_fm_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys_func_multilayer)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),2),...
    'type','spearman');

% MEG
[corr_meg_t_ecm_mst_tdl_resar_subgroup_fm,corr_meg_t_ecm_mst_tdl_resar_subgroup_fm_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_tdl_func_multilayer)), ...
    clin_micro_data(find(select_subs_tdl_for_func_multilayer),10),...
    'type','spearman');
[corr_meg_t_ecm_mst_aprise2040_resar_subgroup_fm,corr_meg_t_ecm_mst_aprise2040_resar_subgroup_fm_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys_func_multilayer)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),12),...
    'type','spearman');
[corr_meg_t_ecm_mst_aprise1st_resar_subgroup_fm,corr_meg_t_ecm_mst_aprise1st_resar_subgroup_fm_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys_func_multilayer)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),2),...
    'type','spearman');

%% Select sample for fMRI / MEG/ multilayer x memory correlations 

%% fMRI 
% RAVLT immediate recall
subs_fmri_15wt_enc = nan(length(data_availability_overview),3);
subs_fmri_15wt_enc(:,1) = data_availability_overview(:,2); 
subs_fmri_15wt_enc(:,2) = data_availability_overview(:,9);
subs_fmri_15wt_enc(:,3) = sum(data_availability_overview(:,[2 9]),2);

select_subs_15wt_enc_for_fmri = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_fmri_15wt_enc(i,2) == 1 && subs_fmri_15wt_enc(i,3) == 2
        select_subs_15wt_enc_for_fmri(i,:) = 1;
    end
end

temp_select_subs_fmri_for_15wt_enc = find(data_availability_overview(:,2));
select_subs_fmri_for_15wt_enc = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_15wt_enc(temp_select_subs_fmri_for_15wt_enc(i,:),2) == 1
        select_subs_fmri_for_15wt_enc(i,:) = 1;
    end
end

% RAVLT delayed recall
subs_fmri_15wt_delrec = nan(length(data_availability_overview),3);
subs_fmri_15wt_delrec(:,1) = data_availability_overview(:,2);
subs_fmri_15wt_delrec(:,2) = data_availability_overview(:,10);
subs_fmri_15wt_delrec(:,3) = sum(data_availability_overview(:,[2 10]),2);

select_subs_15wt_delrec_for_fmri = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_fmri_15wt_delrec(i,2) == 1 && subs_fmri_15wt_delrec(i,3) == 2
        select_subs_15wt_delrec_for_fmri(i,:) = 1;
    end
end

temp_select_subs_fmri_for_15wt_delrec = find(data_availability_overview(:,2));
select_subs_fmri_for_15wt_delrec = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_15wt_delrec(temp_select_subs_fmri_for_15wt_delrec(i,:),2) == 1
        select_subs_fmri_for_15wt_delrec(i,:) = 1;
    end
end

% WMS composite

subs_fmri_wms_alg = nan(length(data_availability_overview),3);
subs_fmri_wms_alg(:,1) = data_availability_overview(:,2);
subs_fmri_wms_alg(:,2) = data_availability_overview(:,12);
subs_fmri_wms_alg(:,3) = sum(data_availability_overview(:,[2 12]),2);

select_subs_wms_alg_for_fmri = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_fmri_wms_alg(i,2) == 1 && subs_fmri_wms_alg(i,3) == 2
        select_subs_wms_alg_for_fmri(i,:) = 1;
    end
end

temp_select_subs_fmri_for_wms_alg = find(data_availability_overview(:,2));
select_subs_fmri_for_wms_alg = zeros(nsubs_fmri,1);
for i = 1:nsubs_fmri
    if subs_fmri_wms_alg(temp_select_subs_fmri_for_wms_alg(i,:),2) == 1
        select_subs_fmri_for_wms_alg(i,:) = 1;
    end
end

%% MEG

% RAVLT immediate recall
subs_meg_15wt_enc = nan(length(data_availability_overview),3);
subs_meg_15wt_enc(:,1) = data_availability_overview(:,3);
subs_meg_15wt_enc(:,2) = data_availability_overview(:,9);
subs_meg_15wt_enc(:,3) = sum(data_availability_overview(:,[3 9]),2);

select_subs_15wt_enc_for_meg = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_meg_15wt_enc(i,2) == 1 && subs_meg_15wt_enc(i,3) == 2
        select_subs_15wt_enc_for_meg(i,:) = 1;
    end
end

temp_select_subs_meg_for_15wt_enc = find(data_availability_overview(:,3));
select_subs_meg_for_15wt_enc = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_15wt_enc(temp_select_subs_meg_for_15wt_enc(i,:),2) == 1
        select_subs_meg_for_15wt_enc(i,:) = 1;
    end
end

% RAVLT delayed recall
subs_meg_15wt_delrec = nan(length(data_availability_overview),3);
subs_meg_15wt_delrec(:,1) = data_availability_overview(:,3);
subs_meg_15wt_delrec(:,2) = data_availability_overview(:,10);
subs_meg_15wt_delrec(:,3) = sum(data_availability_overview(:,[3 10]),2);

select_subs_15wt_delrec_for_meg = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_meg_15wt_delrec(i,2) == 1 && subs_meg_15wt_delrec(i,3) == 2
        select_subs_15wt_delrec_for_meg(i,:) = 1;
    end
end

temp_select_subs_meg_for_15wt_delrec = find(data_availability_overview(:,3));
select_subs_meg_for_15wt_delrec = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_15wt_delrec(temp_select_subs_meg_for_15wt_delrec(i,:),2) == 1
        select_subs_meg_for_15wt_delrec(i,:) = 1;
    end
end

% WMS composite
subs_meg_wms_alg = nan(length(data_availability_overview),3);
subs_meg_wms_alg(:,1) = data_availability_overview(:,3);
subs_meg_wms_alg(:,2) = data_availability_overview(:,12);
subs_meg_wms_alg(:,3) = sum(data_availability_overview(:,[3 12]),2);

select_subs_wms_alg_for_meg = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_meg_wms_alg(i,2) == 1 && subs_meg_wms_alg(i,3) == 2
        select_subs_wms_alg_for_meg(i,:) = 1;
    end
end

temp_select_subs_meg_for_wms_alg = find(data_availability_overview(:,3));
select_subs_meg_for_wms_alg = zeros(nsubs_meg,1);
for i = 1:nsubs_meg
    if subs_meg_wms_alg(temp_select_subs_meg_for_wms_alg(i,:),2) == 1
        select_subs_meg_for_wms_alg(i,:) = 1;
    end
end

%% Multilayer

% RAVLT immediate recall
subs_func_multilayer_15wt_enc = nan(length(data_availability_overview),3);
subs_func_multilayer_15wt_enc(:,1) = data_availability_overview(:,11); % multilayer
subs_func_multilayer_15wt_enc(:,2) = data_availability_overview(:,9); % encoding
subs_func_multilayer_15wt_enc(:,3) = sum(data_availability_overview(:,[11 9]),2);

nsubs_func_multilayer = sum(data_availability_overview(:,11));

select_subs_15wt_enc_for_func_multilayer = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_func_multilayer_15wt_enc(i,2) == 1 && subs_func_multilayer_15wt_enc(i,3) == 2
        select_subs_15wt_enc_for_func_multilayer(i,:) = 1;
    end
end

temp_select_subs_func_multilayer_for_15wt_enc = find(data_availability_overview(:,11));
select_subs_func_multilayer_for_15wt_enc = zeros(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if subs_func_multilayer_15wt_enc(temp_select_subs_func_multilayer_for_15wt_enc(i,:),2) == 1
        select_subs_func_multilayer_for_15wt_enc(i,:) = 1;
    end
end

% RAVLT delayed recall
subs_func_multilayer_15wt_delrec = nan(length(data_availability_overview),3);
subs_func_multilayer_15wt_delrec(:,1) = data_availability_overview(:,11); %multilayer
subs_func_multilayer_15wt_delrec(:,2) = data_availability_overview(:,10); % delrec
subs_func_multilayer_15wt_delrec(:,3) = sum(data_availability_overview(:,[11 10]),2);

select_subs_15wt_delrec_for_func_multilayer = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_func_multilayer_15wt_delrec(i,2) == 1 && subs_func_multilayer_15wt_delrec(i,3) == 2
        select_subs_15wt_delrec_for_func_multilayer(i,:) = 1;
    end
end

temp_select_subs_func_multilayer_for_15wt_delrec = find(data_availability_overview(:,11));
select_subs_func_multilayer_for_15wt_delrec = zeros(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if subs_func_multilayer_15wt_delrec(temp_select_subs_func_multilayer_for_15wt_delrec(i,:),2) == 1
        select_subs_func_multilayer_for_15wt_delrec(i,:) = 1;
    end
end

% WMS composite

subs_func_multilayer_wms_alg = nan(length(data_availability_overview),3);
subs_func_multilayer_wms_alg(:,1) = data_availability_overview(:,11);
subs_func_multilayer_wms_alg(:,2) = data_availability_overview(:,12);
subs_func_multilayer_wms_alg(:,3) = sum(data_availability_overview(:,[11 12]),2);

nsubs_func_multilayer = sum(data_availability_overview(:,11));

select_subs_wms_alg_for_func_multilayer = zeros(length(data_availability_overview),1);
for i = 1:length(data_availability_overview)
    if subs_func_multilayer_wms_alg(i,2) == 1 && subs_func_multilayer_wms_alg(i,3) == 2
        select_subs_wms_alg_for_func_multilayer(i,:) = 1;
    end
end

temp_select_subs_func_multilayer_for_wms_alg = find(data_availability_overview(:,11));
select_subs_func_multilayer_for_wms_alg = zeros(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if subs_func_multilayer_wms_alg(temp_select_subs_func_multilayer_for_wms_alg(i,:),2) == 1
        select_subs_func_multilayer_for_wms_alg(i,:) = 1;
    end
end

%% Correlations 

[corr_fmri_ecm_mst_15wt_enc_resar,corr_fmri_ecm_mst_15wt_enc_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_15wt_enc),:),...
    clin_micro_data(find(select_subs_15wt_enc_for_fmri),13),'type','spearman');
[corr_fmri_ecm_mst_15wt_delrec_resar,corr_fmri_ecm_mst_15wt_delrec_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_15wt_delrec),:),...
    clin_micro_data(find(select_subs_15wt_delrec_for_fmri),14),'type','spearman');
[corr_fmri_ecm_mst_wms_verb_resar,corr_fmri_ecm_mst_wms_verb_resar_p] = ...
    corr(aal_fmri_ecm_mst_resar(find(select_subs_fmri_for_wms_alg),:),...
    clin_micro_data(find(select_subs_wms_alg_for_fmri),15),'type','spearman');

[corr_meg_t_ecm_mst_15wt_enc_resar,corr_meg_t_ecm_mst_15wt_enc_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_15wt_enc),:),...
    clin_micro_data(find(select_subs_15wt_enc_for_meg),13),'type','spearman');
[corr_meg_t_ecm_mst_15wt_delrec_resar,corr_meg_t_ecm_mst_15wt_delrec_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_15wt_delrec),:),...
    clin_micro_data(find(select_subs_15wt_delrec_for_meg),14),'type','spearman');
[corr_meg_t_ecm_mst_wms_verb_resar,corr_meg_t_ecm_mst_wms_verb_resar_p] = ...
    corr(aal_meg_t_ecm_mst_resar(find(select_subs_meg_for_wms_alg),:),...
    clin_micro_data(find(select_subs_wms_alg_for_meg),15),'type','spearman');

[corr_func_multilayer_ecm_mst_15wt_enc_resar,corr_func_multilayer_ecm_mst_15wt_enc_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_15wt_enc),:),...
    clin_micro_data(find(select_subs_15wt_enc_for_func_multilayer),13),'type','spearman');
[corr_func_multilayer_ecm_mst_15wt_delrec_resar,corr_func_multilayer_ecm_mst_15wt_delrec_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_15wt_delrec),:),...
    clin_micro_data(find(select_subs_func_multilayer_for_15wt_delrec),14),'type','spearman');
[corr_func_multilayer_ecm_wms_verb_resar,corr_func_multilayer_ecm_wms_verb_resar_p] = ...
    corr(aal_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_wms_alg),:),...
    clin_micro_data(find(select_subs_wms_alg_for_func_multilayer),15),'type','spearman');


%%
%%%%% ADDITIONAL ANALYSES
%% BNA instead of AAL atlas

%% calculate MST

rois_bna_fmri = 188;

bna_fmri_adjmats_mst = zeros(rois_bna_fmri,rois_bna_fmri,nsubs_fmri);
for i = 1:nsubs_fmri
    temp = bna_fmri_adjmats(:,:,i);
    length(isinf(temp));
    temp(isinf(temp)) = 0;
    [tempmat,~,~] = kruskal_algorithm(temp);
    bna_fmri_adjmats_mst(:,:,i) = tempmat;
end

clear temp*

%% determine network chars

bna_fmri_ecm_mst = nan(nsubs_fmri,rois_bna_fmri);
bna_fmri_k_mst = nan(nsubs_fmri,rois_bna_fmri);
bna_fmri_cc = nan(nsubs_fmri,rois_bna_fmri);
for i = 1:nsubs_fmri
    bna_fmri_ecm_mst(i,:) = (eigenvector_centrality_und(bna_fmri_adjmats_mst(:,:,i)))';
    bna_fmri_k_mst(i,:) = (degrees_und(bna_fmri_adjmats_mst(:,:,i)))';
    bna_fmri_cc(i,:) = (efficiency_wei(bna_fmri_adjmats(:,:,i),2))';
end

%% extract resar
% Left resected area = 74, Right = 75

bna_regions_fmri = lat_fmri;
bna_regions_fmri(lat_fmri == 1,:) = 74;
bna_regions_fmri(lat_fmri == 2,:) = 75;

bna_fmri_ecm_mst_resar = nan(nsubs_fmri,1);
bna_fmri_k_mst_resar = nan(nsubs_fmri,1);
bna_fmri_cc_resar = nan(nsubs_fmri,1);
for i = 1:nsubs_fmri
    bna_fmri_ecm_mst_resar(i,:) = bna_fmri_ecm_mst(i,bna_regions_fmri(i,:));
    bna_fmri_k_mst_resar(i,:) = bna_fmri_k_mst(i,bna_regions_fmri(i,:));
    bna_fmri_cc_resar(i,:) = bna_fmri_cc(i,bna_regions_fmri(i,:));
end

%% correlate with micro

[corr_bna_fmri_cc_tdl_resar,corr_bna_fmri_cc_tdl_resar_p] = ...
    corr(bna_fmri_cc_resar(find(select_subs_fmri_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_fmri),10),'type','spearman');
[corr_bna_fmri_ecm_mst_tdl_resar,corr_bna_fmri_ecm_mst_tdl_resar_p] = ...
    corr(bna_fmri_ecm_mst_resar(find(select_subs_fmri_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_fmri),10),'type','spearman');
[corr_bna_fmri_k_mst_tdl_resar,corr_bna_fmri_k_mst_tdl_resar_p] = ...
    corr(bna_fmri_k_mst_resar(find(select_subs_fmri_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_fmri),10),'type','spearman');

[corr_bna_fmri_cc_aprise1st_resar,corr_bna_fmri_cc_aprise1st_resar_p] = ...
    corr(bna_fmri_cc_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),2),'type','spearman');
[corr_bna_fmri_ecm_mst_aprise1st_resar,corr_bna_fmri_ecm_mst_aprise1st_resar_p] = ...
    corr(bna_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),2),'type','spearman');
[corr_bna_fmri_k_mst_aprise1st_resar,corr_bna_fmri_k_mst_aprise1st_resar_p] = ...
    corr(bna_fmri_k_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),2),'type','spearman');

[corr_bna_fmri_cc_aprise2040_resar,corr_bna_fmri_cc_aprise2040_resar_p] = ...
    corr(bna_fmri_cc_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),12),'type','spearman');
[corr_bna_fmri_ecm_mst_aprise2040_resar,corr_bna_fmri_ecm_mst_aprise2040_resar_p] = ...
    corr(bna_fmri_ecm_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),12),'type','spearman');
[corr_bna_fmri_k_mst_aprise2040_resar,corr_bna_fmri_k_mst_aprise2040_resar_p] = ...
    corr(bna_fmri_k_mst_resar(find(select_subs_fmri_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_fmri),12),'type','spearman');

%% MEG 

rois_bna_meg = 210;

bna_meg_adjmats_t_mst = zeros(rois_bna_meg,rois_bna_meg,nsubs_meg);
for i = 1:nsubs_meg
    temp = bna_meg_adjmats_t(:,:,i);
    length(isinf(temp));
    temp(isinf(temp)) = 0;
    [tempmat,~,~] = kruskal_algorithm(temp);
    bna_meg_adjmats_t_mst(:,:,i) = tempmat;
end

clear temp*

%% determine network chars

bna_meg_t_ecm_mst = nan(nsubs_meg,rois_bna_meg);
bna_meg_t_k_mst = nan(nsubs_meg,rois_bna_meg);
bna_meg_t_cc = nan(nsubs_meg,rois_bna_meg);
for i = 1:nsubs_meg
    bna_meg_t_ecm_mst(i,:) = (eigenvector_centrality_und(bna_meg_adjmats_t_mst(:,:,i)))';
    bna_meg_t_k_mst(i,:) = (degrees_und(bna_meg_adjmats_t_mst(:,:,i)))';
    bna_meg_t_cc(i,:) = (efficiency_wei(bna_meg_adjmats_t(:,:,i),2))';
end

%% extract resar

bna_regions_meg = lat_meg;
bna_regions_meg(lat_meg == 1,:) = 81;
bna_regions_meg(lat_meg == 2,:) = 82;

bna_meg_t_ecm_mst_resar = nan(nsubs_meg,1);
bna_meg_t_k_mst_resar = nan(nsubs_meg,1);
bna_meg_t_cc_resar = nan(nsubs_meg,1);
for i = 1:nsubs_meg
    bna_meg_t_ecm_mst_resar(i,:) = bna_meg_t_ecm_mst(i,bna_regions_meg(i,:));
    bna_meg_t_k_mst_resar(i,:) = bna_meg_t_k_mst(i,bna_regions_meg(i,:));
    bna_meg_t_cc_resar(i,:) = bna_meg_t_cc(i,bna_regions_meg(i,:));
end

%% do correlations

[corr_bna_meg_t_cc_tdl_resar,corr_bna_meg_t_cc_tdl_resar_p] = ...
    corr(bna_meg_t_cc_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');
[corr_bna_meg_t_k_mst_tdl_resar,corr_bna_meg_t_k_mst_tdl_resar_p] = ...
    corr(bna_meg_t_k_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');
[corr_bna_meg_t_ecm_mst_tdl_resar,corr_bna_meg_t_ecm_mst_tdl_resar_p] = ...
    corr(bna_meg_t_ecm_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');

[corr_bna_meg_t_cc_aprise1st_resar,corr_bna_meg_t_cc_aprise1st_resar_p] = ...
    corr(bna_meg_t_cc_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_bna_meg_t_k_mst_aprise1st_resar,corr_bna_meg_t_k_mst_aprise1st_resar_p] = ...
    corr(bna_meg_t_k_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_bna_meg_t_ecm_mst_aprise1st_resar,corr_bna_meg_t_ecm_mst_aprise1st_resar_p] = ...
    corr(bna_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');

[corr_bna_meg_t_cc_aprise2040_resar,corr_bna_meg_t_cc_aprise2040_resar_p] = ...
    corr(bna_meg_t_cc_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');
[corr_bna_meg_t_k_mst_aprise2040_resar,corr_bna_meg_t_k_mst_aprise2040_resar_p] = ...
    corr(bna_meg_t_k_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');
[corr_bna_meg_t_ecm_mst_aprise2040_resar,corr_bna_meg_t_ecm_mst_aprise2040_resar_p] = ...
    corr(bna_meg_t_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');

%% add MEG degree analysis

aal_meg_t_k_mst = nan(nsubs_meg,78);
for i = 1:nsubs_meg
    aal_meg_t_k_mst(i,:) = (degrees_und(aal_meg_adjmats_t_mst(:,:,i)))';
end

aal_meg_t_k_mst_resar = nan(nsubs_meg,1);
for i = 1:nsubs_meg
    aal_meg_t_k_mst_resar(i,:) = aal_meg_t_k_mst(i,aal_regions_meg(i,:));
end

[corr_meg_t_k_mst_aprise1st_resar,corr_meg_t_k_mst_aprise1st_p] = ...
    corr(aal_meg_t_k_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_meg_t_k_mst_aprise2040_resar,corr_meg_t_k_mst_aprise2040_resar_p] = ...
    corr(aal_meg_t_k_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');
[corr_meg_t_k_mst_tdl_resar,corr_meg_t_k_mst_tdl_resar_p] = ...
    corr(aal_meg_t_k_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');

%% also do ecm mst analysis for a1 and bb bands

aal_meg_a1_ecm_mst = nan(nsubs_meg,78);
aal_meg_bb_ecm_mst = nan(nsubs_meg,78);
for i = 1:nsubs_meg
    aal_meg_a1_ecm_mst(i,:) = (eigenvector_centrality_und(aal_meg_adjmats_a1_mst(:,:,i)))';
    aal_meg_bb_ecm_mst(i,:) = (eigenvector_centrality_und(aal_meg_adjmats_bb_mst(:,:,i)))';
end

aal_meg_a1_ecm_mst_resar = nan(nsubs_meg,1);
aal_meg_bb_ecm_mst_resar = nan(nsubs_meg,1);
for i = 1:nsubs_meg
    aal_meg_a1_ecm_mst_resar(i,:) = aal_meg_a1_ecm_mst(i,aal_regions_meg(i,:));
    aal_meg_bb_ecm_mst_resar(i,:) = aal_meg_bb_ecm_mst(i,aal_regions_meg(i,:));
end

[corr_meg_a1_ecm_mst_aprise1st_resar,corr_meg_a1_ecm_mst_aprise1st_resar_p] = ...
    corr(aal_meg_a1_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_meg_a1_ecm_mst_aprise2040_resar,corr_meg_a1_ecm_mst_aprise2040_resar_p] = ...
    corr(aal_meg_a1_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');
[corr_meg_a1_ecm_mst_tdl_resar,corr_meg_a1_ecm_mst_tdl_resar_p] = ...
    corr(aal_meg_a1_ecm_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');

[corr_meg_bb_ecm_mst_aprise1st_resar,corr_meg_bb_ecm_mst_aprise1st_p] = ...
    corr(aal_meg_bb_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),2),'type','spearman');
[corr_meg_bb_ecm_mst_aprise2040_resar,corr_meg_bb_ecm_mst_aprise2040_resar_p] = ...
    corr(aal_meg_bb_ecm_mst_resar(find(select_subs_meg_for_ephys),:),...
    clin_micro_data(find(select_subs_ephys_for_meg),12),'type','spearman');
[corr_meg_bb_ecm_mst_tdl_resar,corr_meg_bb_ecm_mst_tdl_resar_p] = ...
    corr(aal_meg_bb_ecm_mst_resar(find(select_subs_meg_for_tdl),:),...
    clin_micro_data(find(select_subs_tdl_for_meg),10),'type','spearman');

%% BNA multilayer

bna_regions_func_multilayer = nan(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    if lat_func_multilayer(i,:) == 1 % 1 = left
        bna_regions_func_multilayer(i,:) = 74;
    elseif lat_func_multilayer(i,:) == 2 % 2 = right
        bna_regions_func_multilayer(i,:) = 75;
    end
end

%% extract bc of resected area

bna_func_multilayer_ecm_resar = nan(nsubs_func_multilayer,1);
for i = 1:nsubs_func_multilayer
    bna_func_multilayer_ecm_resar(i,:) = bna_func_multilayer_ecm(bna_regions_func_multilayer(i,:),i);
end

%% correlate this to micro

[corr_bna_multilayer_ecm_mst_aprise1st_resar,corr_bna_multilayer_ecm_mst_aprise1st_resar_p] = ...
    corr(bna_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),2),...
    'type','spearman');
[corr_bna_multilayer_ecm_mst_aprise2040_resar,corr_bna_multilayer_ecm_mst_aprise2040_resar_p] = ...
    corr(bna_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_ephys)), ...
    clin_micro_data(find(select_subs_ephys_for_func_multilayer),12),...
    'type','spearman');
[corr_bna_multilayer_ecm_tdl_resar,corr_bna_multilayer_ecm_tdl_resar_p] = ...
    corr(bna_func_multilayer_ecm_resar(find(select_subs_func_multilayer_for_tdl)), ...
    clin_micro_data(find(select_subs_tdl_for_func_multilayer),10),...
    'type','spearman');


