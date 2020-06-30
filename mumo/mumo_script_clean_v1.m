%% mumo_script_clean_v1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 14.05.2020
%%% 01.06.2020
%%% L BREEDT
%
% This script contains all of the code that was used in the MuMoBrain
% project.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PART 1: CALCULATING CONNECTIVITY
addpath('/mnt/anw-gold/MUMO/6_projects/breedt_mumo_pilot/b_analyses/scripts/');
addpath('/mnt/anw-gold/MUMO/5_proc_data/');


% -------------------------------- MEG ------------------------------------

% load cases
cd('/mnt/anw-gold/MUMO/5_proc_data/meg_BNA_asciis_version1');
fname = dir('mumo*'); % all directories of included subjects

% set variables
Fs = 1250; % sampling frequency
%window_size = 4096;
%overlap = window_size/2;
%num_ffts = window_size;

epoch = 1:22; % number of epochs to analyze
epoch_length = 1:4096:16384; % split epochs?
nrois = 210; % number of rois to analyze (BNA = 246, cortical BNA = 210)

% pre-allocate matrices
pli_delta_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois);
pli_theta_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_alpha1_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_alpha2_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois);
pli_beta_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_gamma_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), nrois, nrois);
%PLI_HC_broad_epoch = zeros(numel(fname), numel(epoch), numel(epoch_length), 78, 78);

% loop over all cases and epochs and calculate PLI
tic
for k = 1:length(fname) % k = case #
    cd(['/mnt/anw-gold/MUMO/5_proc_data/meg_BNA_asciis_version1/' fname(k).name '/OD1/']) % first measurement only
    files = dir('*.asc'); % all epochs of 1 subject
    fprintf(1, 'Now calculating PLI for sub %s\n', num2str(k))
    for i = 1:length(epoch)       % i = epoch #
        data = importdata(fullfile(pwd, files(i).name)); % import epoch data
        data = data(:,1:nrois); % subset of specified rois data only
        for j = 1:numel(epoch_length)
            data_epoch = data(epoch_length(j):epoch_length(j)+4095,:);
            
            pli_delta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, .5, 4, nrois));
            pli_theta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 4, 8, nrois)); 
            pli_alpha1_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 8, 10, nrois));
            pli_alpha2_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 10, 13, nrois));
            pli_beta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 13, 30, nrois));
            pli_gamma_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 30, 48, nrois));
            %PLI_HC_broad_epoch (k,i,j,:,:) = pli_matteo(fft_filt(data_epoch, Fs, .5, 48));         
        end
    end
end
toc

% average over epochs (and get 3D array, SUBxROIxROI)
x=squeeze(mean(pli_delta_epoch,3));
pli_delta_full_raw=squeeze(mean(x,2));

x=squeeze(mean(pli_theta_epoch,3));
pli_theta_full_raw=squeeze(mean(x,2));

x=squeeze(mean(pli_alpha1_epoch,3));
pli_alpha1_full_raw=squeeze(mean(x,2));

x=squeeze(mean(pli_alpha2_epoch,3));
pli_alpha2_full_raw=squeeze(mean(x,2));

x=squeeze(mean(pli_beta_epoch,3));
pli_beta_full_raw=squeeze(mean(x,2));

x=squeeze(mean(pli_gamma_epoch,3));
pli_gamma_full_raw=squeeze(mean(x,2));

% save raw pli matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/full_raw.mat','*full_raw')


% ------------------------------- fMRI ------------------------------------

% specify filepath & set variables
path = '/mnt/anw-gold/MUMO/5_proc_data/fmri_BNA_timeseries_version1/';
filepat = fullfile(path, 'mumo*');

subs = dir(filepat);
nsubs = length(subs);
nrois = 225;

% loop over all cases and calculate pearson's correlations
tic
fmri_full_raw = zeros(nsubs, nrois, nrois);
for k = 1:nsubs
    basefname = subs(k).name;
    fullfname = fullfile(path, basefname);
    fprintf(1, 'Currently calculating pearson correlations for %s\n', basefname)
    sub = dlmread(fullfname);    
    adjmat = abs(corr(sub));
    adjmat(logical(eye(size(adjmat)))) = 0;
    adjmat(isnan(adjmat)) = 0;
    
    fmri_full_raw(k,:,:) = adjmat;
end
toc

% save raw fmri matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/full_raw.mat','fmri_full_raw','-append')

% remove subcortical rois
fmri_full_raw(:,211:225,:) = [];
fmri_full_raw(:,:,211:225) = [];


% -------------------------------- DWI ------------------------------------

% specify filepath & set variables
path = '/mnt/anw-gold/MUMO/5_proc_data/dwi_BNA/';
filepat = fullfile(path, '*BNA.csv');

subs = dir(filepat); 
nsubs = length(subs);
nrois = 224; 

dwi_full_raw = zeros(nsubs, nrois, nrois);
for k = 1:nsubs
    basefname = subs(k).name;
    fullfname = fullfile(path, basefname);
    
    sub = dlmread(fullfname);
    
    adjmat = triu(sub) + triu(sub)';
    adjmat = adjmat.*~eye(size(adjmat));
    
    dwi_full_raw(k,:,:) = adjmat;
end

% save raw dwi matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/full_raw.mat','dwi_full_raw','-append')

% remove subcortical rois
dwi_full_raw(:,211:224,:) = [];
dwi_full_raw(:,:,211:224) = [];


% ---------------------------- remove ROIs --------------------------------

% find + remove fmri + dwi empty regions
[fmri_full_clean, deleted_regions_fmri, dwi_full_clean, deleted_regions_dwi] = remove_empty_regions(fmri_full_raw, dwi_full_raw);

% concatenate fmri + dwi deleted regions
deleted_regions = vertcat(deleted_regions_fmri, deleted_regions_dwi);
deleted_regions = sort(deleted_regions);

% remove empty regions from MEG
pli_delta_full_clean(:,:,deleted_regions)=[];
pli_delta_full_clean(:,deleted_regions,:)=[];

pli_theta_full_clean(:,:,deleted_regions)=[];
pli_theta_full_clean(:,deleted_regions,:)=[];

pli_alpha1_full_clean(:,:,deleted_regions)=[];
pli_alpha1_full_clean(:,deleted_regions,:)=[];

pli_alpha2_full_clean(:,:,deleted_regions)=[];
pli_alpha2_full_clean(:,deleted_regions,:)=[];

pli_beta_full_clean(:,:,deleted_regions)=[];
pli_beta_full_clean(:,deleted_regions,:)=[];

pli_gamma_full_clean(:,:,deleted_regions)=[];
pli_gamma_full_clean(:,deleted_regions,:)=[];

% save clean matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/full_clean.mat','*full_clean')


% ----------------------------- normalize ---------------------------------

% MEG
pli_delta_full_norm = zeros(nsubs,nrois,nrois);
pli_theta_full_norm = zeros(nsubs,nrois,nrois);
pli_alpha1_full_norm = zeros(nsubs,nrois,nrois);
pli_alpha2_full_norm = zeros(nsubs,nrois,nrois);
pli_beta_full_norm = zeros(nsubs,nrois,nrois);
pli_gamma_full_norm = zeros(nsubs,nrois,nrois);
for sub = 1:nsubs
    M = squeeze(pli_delta_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_delta_full_norm(sub,:,:) = weight_conversion(M,'normalize');
    
    M = squeeze(pli_theta_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_theta_full_norm(sub,:,:) = weight_conversion(M,'normalize');
    
    M = squeeze(pli_alpha1_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_alpha1_full_norm(sub,:,:) = weight_conversion(M,'normalize');
    
    M = squeeze(pli_alpha2_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_alpha2_full_norm(sub,:,:) = weight_conversion(M,'normalize');
    
    M = squeeze(pli_beta_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_beta_full_norm(sub,:,:) = weight_conversion(M,'normalize');
    
    M = squeeze(pli_gamma_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    pli_gamma_full_norm(sub,:,:) = weight_conversion(M,'normalize');
end

% fMRI
fmri_full_norm = zeros(nsubs,nrois,nrois);
for sub = 1:nsubs
    M = squeeze(fmri_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    fmri_full_norm(sub,:,:) = weight_conversion(M,'normalize');
end

% DWI
dwi_full_norm = zeros(nsubs,nrois,nrois);
for sub = 1:nsubs
    M = squeeze(dwi_full_clean(sub,:,:));
    M = weight_conversion(M,'autofix');
    dwi_full_norm(sub,:,:) = weight_conversion(M,'normalize');
end

% save normalized matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/clean_normalized.mat','*full_norm')


% -------------------------------- MST ------------------------------------

% set variables
nsubs = 33;
nrois = 197;

% MEG
% pre-allocate matrices
mst_pli_del = zeros(nsubs, nrois, nrois);
mst_pli_the = zeros(nsubs, nrois, nrois);
mst_pli_al1 = zeros(nsubs, nrois, nrois);
mst_pli_al2 = zeros(nsubs, nrois, nrois);
mst_pli_bet = zeros(nsubs, nrois, nrois);
mst_pli_gam = zeros(nsubs, nrois, nrois);
% construct mst
for sub = 1:nsubs
    mst_pli_del(sub,:,:) = kruskal_algorithm(squeeze(pli_delta_full_raw(sub,:,:)));
    mst_pli_the(sub,:,:) = kruskal_algorithm(squeeze(pli_theta_full_raw(sub,:,:)));
    mst_pli_al1(sub,:,:) = kruskal_algorithm(squeeze(pli_alpha1_full_raw(sub,:,:)));
    mst_pli_al2(sub,:,:) = kruskal_algorithm(squeeze(pli_alpha2_full_raw(sub,:,:)));
    mst_pli_bet(sub,:,:) = kruskal_algorithm(squeeze(pli_beta_full_raw(sub,:,:)));
    mst_pli_gam(sub,:,:) = kruskal_algorithm(squeeze(pli_gamma_full_raw(sub,:,:)));
end

% fMRI
% pre-allocate matrices
mst_fmri = zeros(nsubs, nrois, nrois);
% construct mst
for sub = 1:nsubs
    mst_fmri(sub,:,:) = kruskal_algorithm(squeeze(fmri_full_raw(sub,:,:)));
end

% DWI
% pre-allocate matrices
mst_dwi = zeros(nsubs, nrois, nrois);
% construct mst
for sub = 1:nsubs
    mst_dwi(sub,:,:) = kruskal_algorithm(squeeze(dwi_full_raw(sub,:,:)));
end

% save MSTs
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/clean_mst.mat','mst*')



%% PART 2: RANDOMIZING DATA

% ---------------------------- weighted -----------------------------------

% set variables
nsubs = 33;
nrois = 197;

% pre-allocate matrices
rand_fmri = zeros(nsubs,nrois,nrois);

rand_pli_del = zeros(nsubs,nrois,nrois);
rand_pli_the = zeros(nsubs,nrois,nrois);
rand_pli_al1 = zeros(nsubs,nrois,nrois);
rand_pli_al2 = zeros(nsubs,nrois,nrois);
rand_pli_bet = zeros(nsubs,nrois,nrois);
rand_pli_gam = zeros(nsubs,nrois,nrois);

rand_dwi = zeros(nsubs,nrois,nrois);

% loop & shuffle
tic
for sub = 1:nsubs
    fprintf(1, 'Now shuffling sub %s!\n', num2str(sub))
    rand_fmri(sub,:,:) = null_model_und_sign(squeeze(fmri_full_norm(sub,:,:)),1000,1);
    
    rand_pli_del(sub,:,:) = null_model_und_sign(squeeze(pli_delta_full_raw(sub,:,:)),1000,1);
    rand_pli_the(sub,:,:) = null_model_und_sign(squeeze(pli_theta_full_raw(sub,:,:)),1000,1);
    rand_pli_al1(sub,:,:) = null_model_und_sign(squeeze(pli_alpha1_full_raw(sub,:,:)),1000,1);
    rand_pli_al2(sub,:,:) = null_model_und_sign(squeeze(pli_alpha2_full_raw(sub,:,:)),1000,1);
    rand_pli_bet(sub,:,:) = null_model_und_sign(squeeze(pli_beta_full_raw(sub,:,:)),1000,1);
    rand_pli_gam(sub,:,:) = null_model_und_sign(squeeze(pli_gamma_full_raw(sub,:,:)),1000,1);
    
    rand_dwi(sub,:,:) = null_model_und_sign(squeeze(dwi_full_norm(sub,:,:)),1000,1);
end
toc

% save randomized matrices
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/clean_randomized.mat','rand*')


% -------------------------------- MST ------------------------------------

% pre-allocate matrices
rand_mst_fmri = zeros(nsubs,nrois,nrois);

rand_mst_pli_del = zeros(nsubs,nrois,nrois);
rand_mst_pli_the = zeros(nsubs,nrois,nrois);
rand_mst_pli_al1 = zeros(nsubs,nrois,nrois);
rand_mst_pli_al2 = zeros(nsubs,nrois,nrois);
rand_mst_pli_bet = zeros(nsubs,nrois,nrois);
rand_mst_pli_gam = zeros(nsubs,nrois,nrois);

rand_mst_dwi = zeros(nsubs,nrois,nrois);

% construct MST from random matrices
for sub = 1:nsubs
    fprintf(1, 'Now calculating MST for sub %s!\n', num2str(sub))
    rand_mst_fmri(sub,:,:) = kruskal_algorithm(squeeze(rand_fmri(sub,:,:)));
    
    rand_mst_pli_del(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_del(sub,:,:)));
    rand_mst_pli_the(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_the(sub,:,:)));
    rand_mst_pli_al1(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_al1(sub,:,:)));
    rand_mst_pli_al2(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_al2(sub,:,:)));
    rand_mst_pli_bet(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_bet(sub,:,:)));
    rand_mst_pli_gam(sub,:,:) = kruskal_algorithm(squeeze(rand_pli_gam(sub,:,:)));
    
    rand_mst_dwi(sub,:,:) = kruskal_algorithm(squeeze(rand_dwi(sub,:,:)));
end

% save randomized MSTs
%save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/clean_randomized_mst.mat','rand_mst*')



%% PART 3: SUPRA-ADJACENCY MATRICES

% set variables
nlrs = 8;
id = nrois:nrois:nrois*nlrs;

% ---------------------------- weighted -----------------------------------

% pre-allocate supra-adjacency matrices
supra_weighted_full = zeros(nrois*nlrs, nrois*nlrs, nsubs);
rand_supra_weighted = zeros(nrois*nlrs, nrois*nlrs, nsubs);

% construct weighted supra-adjacency matrix
for c = 1:nsubs
    fprintf(1, 'Now constructing supra_weighted for sub %s!\n', num2str(c))
    
    supratmp = blkdiag(squeeze(fmri_full_norm(c,:,:)), squeeze(pli_delta_full_norm(c,:,:)), squeeze(pli_theta_full_norm(c,:,:)), squeeze(pli_alpha1_full_norm(c,:,:)), squeeze(pli_alpha2_full_norm(c,:,:)), squeeze(pli_beta_full_norm(c,:,:)), squeeze(pli_gamma_full_norm(c,:,:)), squeeze(dwi_full_norm(c,:,:)));
    for i = 1:length(id)
        supratmp(id(i)*size(supratmp,1)+1:size(supratmp,1)+1:end) = 1;
        supratmp(id(i)+1:size(supratmp, 1)+1:1+size(supratmp, 1)*min(size(supratmp, 1)-id(i),size(supratmp, 2))) = 1;
    end
    supra_weighted_full(:,:,c) = supratmp;
end

% save weighted supra-adjacency matrix
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/supra_weighted_full.mat','supra_weighted_full')

% construct randomized weighted supra-adjacency matrix
for c = 1:nsubs
    fprintf(1, 'Now constructing rand_supra_weighted for sub %s!\n', num2str(c))
    
    supratmp = blkdiag(squeeze(rand_fmri(c,:,:)), squeeze(rand_pli_del(c,:,:)), squeeze(rand_pli_the(c,:,:)), squeeze(rand_pli_al1(c,:,:)), squeeze(rand_pli_al2(c,:,:)), squeeze(rand_pli_bet(c,:,:)), squeeze(rand_pli_gam(c,:,:)), squeeze(rand_dwi(c,:,:)));
    for i = 1:length(id)
        supratmp(id(i)*size(supratmp,1)+1:size(supratmp,1)+1:end) = 1;
        supratmp(id(i)+1:size(supratmp, 1)+1:1+size(supratmp, 1)*min(size(supratmp, 1)-id(i),size(supratmp, 2))) = 1;
    end
    rand_supra_weighted(:,:,c) = supratmp;
end

% save randomized weighted supra-adjacency matrix
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/rand_supra_weighted.mat','rand_supra_weighted')


% -------------------------------- MST ------------------------------------

% pre-allocate supra-adjacency matrices
supra_mst_full = zeros(nrois*nlrs, nrois*nlrs, nsubs);
rand_supra_mst = zeros(nrois*nlrs, nrois*nlrs, nsubs);

% construct mst supra-adjacency matrix
for c = 1:nsubs
    fprintf(1, 'Now constructing supra_mst for sub %s!\n', num2str(c))
    
    supratmp = blkdiag(squeeze(mst_fmri(c,:,:)), squeeze(mst_pli_del(c,:,:)), squeeze(mst_pli_the(c,:,:)), squeeze(mst_pli_al1(c,:,:)), squeeze(mst_pli_al2(c,:,:)), squeeze(mst_pli_bet(c,:,:)), squeeze(mst_pli_gam(c,:,:)), squeeze(mst_dwi(c,:,:)));
    for i = 1:length(id)
        supratmp(id(i)*size(supratmp,1)+1:size(supratmp,1)+1:end) = 1;
        supratmp(id(i)+1:size(supratmp, 1)+1:1+size(supratmp, 1)*min(size(supratmp, 1)-id(i),size(supratmp, 2))) = 1;
    end
    supra_mst_full(:,:,c) = supratmp;
end

% save mst supra-adjacency matrix
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/supra_mst_full.mat','supra_mst_full')

% construct randomized mst supra-adjacency matrix
for c = 1:nsubs
    fprintf(1, 'Now constructing rand_supra_mst for sub %s!\n', num2str(c))
    
    supratmp = blkdiag(squeeze(rand_mst_fmri(c,:,:)), squeeze(rand_mst_pli_del(c,:,:)), squeeze(rand_mst_pli_the(c,:,:)), squeeze(rand_mst_pli_al1(c,:,:)), squeeze(rand_mst_pli_al2(c,:,:)), squeeze(rand_mst_pli_bet(c,:,:)), squeeze(rand_mst_pli_gam(c,:,:)), squeeze(rand_mst_dwi(c,:,:)));
    for i = 1:length(id)
        supratmp(id(i)*size(supratmp,1)+1:size(supratmp,1)+1:end) = 1;
        supratmp(id(i)+1:size(supratmp, 1)+1:1+size(supratmp, 1)*min(size(supratmp, 1)-id(i),size(supratmp, 2))) = 1;
    end
    rand_supra_mst(:,:,c) = supratmp;
end

% save randomized mst supra-adjacency matrix
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/rand_supra_mst.mat','rand_supra_mst')



%% PART 4: NETWORK MEASURES

% set variables
N = 197;

% determine indices of areas belonging to FPN
BNA = 1:210;
BNA(:,deleted_regions) = []; % deleted_regions contains the regions removed in PART 1
FPN =  [17 18 19 20 21 22 29 30 31 32 99 100 137 138 147 148 177 178]; % actual FPN regions
bool = ismember(BNA, FPN);
fpn_id = find(bool); % mumo FPN regions 


% ---------------------------- monolayer ----------------------------------

% pre-allocate matrices
ec_meg_delta_raw = zeros(33, 197);
ec_meg_theta_raw = zeros(33, 197);
ec_meg_alpha1_raw = zeros(33, 197);
ec_meg_alpha2_raw = zeros(33, 197);
ec_meg_beta_raw = zeros(33, 197);
ec_meg_gamma_raw = zeros(33, 197);

ec_meg_delta_fpn = zeros(33, 18);
ec_meg_theta_fpn = zeros(33, 18);
ec_meg_alpha1_fpn = zeros(33, 18);
ec_meg_alpha2_fpn = zeros(33, 18);
ec_meg_beta_fpn = zeros(33, 18);
ec_meg_gamma_fpn = zeros(33, 18);

ec_fmri_raw = zeros(33, 197);
ec_fmri_fpn = zeros(33, 18);

ec_dwi_raw = zeros(33, 197);
ec_dwi_fpn = zeros(33, 18);

% caclulate nodal EC
for sub = 1:33
    ec_meg_delta_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_del(sub,:,:)));
    ec_meg_theta_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_the(sub,:,:)));
    ec_meg_alpha1_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_al1(sub,:,:)));
    ec_meg_alpha2_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_al2(sub,:,:)));
    ec_meg_beta_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_bet(sub,:,:)));
    ec_meg_gamma_raw(sub,:) = eigenvector_centrality_und(squeeze(mst_pli_gam(sub,:,:)));
    
    ec_fmri_raw(sub,:) = eigenvector_centrality_und(squeeze(fmri_full_norm(sub,:,:)));
    
    ec_dwi_raw(sub,:) = eigenvector_centrality_und(squeeze(dwi_norm(sub,:,:)));
end

% normalize EC by dividing by the max to conform with multilayer scaling
ec_meg_delta_norm = zeros(33, 197);
ec_meg_theta_norm = zeros(33, 197);
ec_meg_alpha1_norm = zeros(33, 197);
ec_meg_alpha2_norm = zeros(33, 197);
ec_meg_beta_norm = zeros(33, 197);
ec_meg_gamma_norm = zeros(33, 197);

ec_fmri_norm = zeros(33, 197);

ec_dwi_norm = zeros(33, 197);

for sub = 1:33
    ec_meg_delta_norm(sub,:) = squeeze(ec_meg_delta_raw(sub,:)) / max(squeeze(ec_meg_delta_raw(sub,:)));
    ec_meg_theta_norm(sub,:) = squeeze(ec_meg_theta_raw(sub,:)) / max(squeeze(ec_meg_theta_raw(sub,:)));
    ec_meg_alpha1_norm(sub,:) = squeeze(ec_meg_alpha1_raw(sub,:)) / max(squeeze(ec_meg_alpha1_raw(sub,:)));
    ec_meg_alpha2_norm(sub,:) = squeeze(ec_meg_alpha2_raw(sub,:)) / max(squeeze(ec_meg_alpha2_raw(sub,:)));
    ec_meg_beta_norm(sub,:) = squeeze(ec_meg_beta_raw(sub,:)) / max(squeeze(ec_meg_beta_raw(sub,:)));
    ec_meg_gamma_norm(sub,:) = squeeze(ec_meg_gamma_raw(sub,:)) / max(squeeze(ec_meg_gamma_raw(sub,:)));
    
    ec_fmri_norm(sub,:) = squeeze(ec_fmri_raw(sub,:)) / max(squeeze(ec_fmri_raw(sub,:)));
    
    ec_dwi_norm(sub,:) = squeeze(ec_dwi_raw(sub,:)) / max(squeeze(ec_dwi_raw(sub,:)));
end

% save normalized monolayer nodal EC
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/monolayer_nodal_ec.mat','ec_*_norm')

% extract EC of FPN
for sub = 1:33
    tmp = squeeze(ec_meg_delta_norm(sub,:));
    ec_meg_delta_fpn(sub,:) = tmp(fpn_id);
    tmp = squeeze(ec_meg_theta_norm(sub,:));
    ec_meg_theta_fpn(sub,:) = tmp(fpn_id);
    tmp = squeeze(ec_meg_alpha1_norm(sub,:));
    ec_meg_alpha1_fpn(sub,:) = tmp(fpn_id);
    tmp = squeeze(ec_meg_alpha2_norm(sub,:));
    ec_meg_alpha2_fpn(sub,:) = tmp(fpn_id);
    tmp = squeeze(ec_meg_beta_norm(sub,:));
    ec_meg_beta_fpn(sub,:) = tmp(fpn_id);
    tmp = squeeze(ec_meg_gamma_norm(sub,:));
    ec_meg_gamma_fpn(sub,:) = tmp(fpn_id);
    
    tmp = squeeze(ec_fmri_norm(sub,:));
    ec_fmri_fpn(sub,:) = tmp(fpn_id);
    
    tmp = squeeze(ec_dwi_norm(sub,:));
    ec_dwi_fpn(sub,:) = tmp(fpn_id);
end

% calculate average EC of FPN
ec_meg_delta_fpn_mean = squeeze(mean(ec_meg_delta_fpn, 2));
ec_meg_theta_fpn_mean = squeeze(mean(ec_meg_theta_fpn, 2));
ec_meg_alpha1_fpn_mean = squeeze(mean(ec_meg_alpha1_fpn, 2));
ec_meg_alpha2_fpn_mean = squeeze(mean(ec_meg_alpha2_fpn, 2));
ec_meg_beta_fpn_mean = squeeze(mean(ec_meg_beta_fpn, 2));
ec_meg_gamma_fpn_mean = squeeze(mean(ec_meg_gamma_fpn, 2));

ec_fmri_fpn_mean = squeeze(mean(ec_fmri_fpn, 2));

ec_dwi_fpn_mean = squeeze(mean(ec_dwi_fpn, 2));


% ---------------------------- multilayer ---------------------------------

% load multilayer EC
multilayer_ec_raw = readtable('EC_No_Mask_Group_MST_Multi_layer_real_.csv', 'HeaderLines', 1);

% place in SUBSxVALUES matrix
id_ec = 1:197:6502;
multilayer_ec_mat = zeros(33, 197);
for i = 1:33
    multilayer_ec_mat(i,:) = multilayer_ec_raw{(id_ec(i):id_ec(i+1)-1), 2}';
end

% save multilayer nodal EC
save('/data/MUMO/6_projects/breedt_mumo_pilot/b_analyses/output/multilayer_nodal_ec.mat','multilayer_ec_mat')

% extract FPN
multilayer_ec_fpn_mat = zeros(33, 18);
for sub = 1:33
    tmp = squeeze(multilayer_ec_mat(sub,:));
    multilayer_ec_fpn_mat(sub,:) = tmp(fpn_id);
end

% calculate average EC of FPN
multilayer_ec_fpn = squeeze(mean(multilayer_ec_fpn_mat, 2));


% ------------------------------- export -----------------------------------
%%% write data to a table, containing subject number in the first column,
%%% multilayer EC in the 2nd and 3d columns, and the monolayer EC values in
%%% the following columns; then write this table to a .csv file for easy
%%% importation in SPSS

% create vector containing subject names (corresponding to SPSS file)
PPN = [110002, 110003, 110004, 110005, 110006, 110008, 110009, 110010, 110011, 110012, 110013, 110014, 110015, 110016, 110017, 110018, 110020, 110021, 110023, 110024, 110025, 110027, 110028, 110029, 110030, 110031, 110032, 110033, 110034, 110035, 110037, 110038, 110039]';

% create table containing network measures
ec_table = table(PPN, multilayer_ec_fpn, ec_meg_delta_fpn_mean, ec_meg_theta_fpn_mean, ec_meg_alpha1_fpn_mean, ec_meg_alpha2_fpn_mean, ec_meg_beta_fpn_mean, ec_meg_gamma_fpn_mean, ec_fmri_fpn_mean, ec_dwi_fpn_mean);
writetable(ec_table, '/mnt/anw-gold/MUMO/6_projects/breedt_mumo_pilot/e_miscellaneous/ec_vx.csv')



