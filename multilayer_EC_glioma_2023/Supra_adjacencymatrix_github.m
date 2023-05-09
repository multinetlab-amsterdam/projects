%% PART 1: CALCULATING CONNECTIVITY
%indicate where to find your data (scripts + MEG files)
addpath('/mnt/resource/m.vanlingen/m2b/'); 
addpath('/data/KNW/NO-cohorten/Scans/'); 
addpath('/mnt/resource/m.vanlingen/m2b/2019_03_03_BCT'); %this folder contains all functions that are needed

% -------------------------------- MEG ------------------------------------

% load cases -> refer to subject folders with all MEG asciis.
cd('/data/KNW/NO-cohorten/Scans/');
fname = dir('subjects_alltimepoints.txt'); % textfile indicating all directories of included subjects at T1 and T2
addpath('/mnt/resource/m.vanlingen/m2b/'); 
addpath('/data/KNW/NO-cohorten/Scans/');  


% -------------------------------- MEG ------------------------------------

% load cases -> direct to folders with MEG files (asciis) of all subjects

cd('/data/KNW/NO-cohorten/Scans/');
textfile = 'subjects_alltimepoints.txt';
fname = readtable(textfile, 'ReadVariableNames', false); 
fname = table2cell(fname);

name=fname(:,1);
timepoint= fname(:,2);



nsubs =  length(fname); 

% set variables
Fs = 1250; % sampling frequency


epoch = 1:15; % number of epochs to analyze - same for all participants; at least 13 epochs. 
epoch_length = 1:8192:16384; % split epochs samples by 2. 
nrois = 78; % number of rois to analyze (AAL=78 or alternatively BNA = 246, cortical BNA = 210)

% pre-allocate matrices
pli_delta_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_theta_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_alpha1_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_alpha2_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois);
pli_beta_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois); 
pli_gamma_epoch = zeros(numel(name), numel(epoch), numel(epoch_length), nrois, nrois);


% loop over all cases and epochs and calculate PLI to make the adjacency matrix
tic
for k = 1:nsubs 
    cd(['/data/KNW/NO-cohorten/Scans/' name{k} '/meg/' timepoint{k} '/AAL']) % looping over name & timepoint at the same time. 
    files = dir('*.asc'); % all epochs of 1 subject
    fprintf(1, 'Now calculating PLI for sub %s\n', num2str(k))
    for i = 1:length(epoch)       % i = epoch # data import per epoch
        data = importdata(fullfile(pwd, files(i).name)); % import epoch data
        data = data(:,1:nrois); % subset of specified rois data only
        for j = 1:numel(epoch_length) 
            data_epoch = data(epoch_length(j):epoch_length(j)+8191,:);
            
            pli_delta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, .5, 4, nrois)); 
            pli_theta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 4, 8, nrois)); 
            pli_alpha1_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 8, 10, nrois));
            pli_alpha2_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 10, 13, nrois));
            pli_beta_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 13, 30, nrois));
            pli_gamma_epoch (k,i,j,:,:) = pli_matteo(fft_filt_BNA(data_epoch, Fs, 30, 48, nrois));        
        end
    end
end
toc

% average over epochs (and get 3D array, SUBxROIxROI). Output = 1 adjacency matrix per participant.

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

% ----------------------------- normalize-> do autofix---------------------------------


% set variables
nrois = length(pli_delta_full_raw(1,1,:));

% MEG
pli_delta_full_autofix = zeros(nsubs,nrois,nrois);
pli_theta_full_autofix = zeros(nsubs,nrois,nrois);
pli_alpha1_full_autofix = zeros(nsubs,nrois,nrois);
pli_alpha2_full_autofix = zeros(nsubs,nrois,nrois);
pli_beta_full_autofix = zeros(nsubs,nrois,nrois);
pli_gamma_full_autofix = zeros(nsubs,nrois,nrois);

for sub = 1:nsubs
    M = squeeze(pli_delta_full_raw(sub,:,:));
    pli_delta_full_autofix(sub,:,:) = weight_conversion(M,'autofix');
  
    M = squeeze(pli_theta_full_raw(sub,:,:));
    pli_theta_full_autofix(sub,:,:)= weight_conversion(M,'autofix');
    
    M = squeeze(pli_alpha1_full_raw(sub,:,:));
    pli_alpha1_full_autofix(sub,:,:)= weight_conversion(M,'autofix');
    
    M = squeeze(pli_alpha2_full_raw(sub,:,:));
    pli_alpha2_full_autofix(sub,:,:)= weight_conversion(M,'autofix');
    
    M = squeeze(pli_beta_full_raw(sub,:,:));
    pli_beta_full_autofix(sub,:,:)= weight_conversion(M,'autofix'); 
    
    M = squeeze(pli_gamma_full_raw(sub,:,:));
    pli_gamma_full_autofix(sub,:,:)= weight_conversion(M,'autofix');
    
end



% --------------------pre-allocate matrices----------------------------------------------------------------
mst_pli_del = zeros(nsubs, nrois, nrois);
mst_pli_the = zeros(nsubs, nrois, nrois);
mst_pli_al1 = zeros(nsubs, nrois, nrois);
mst_pli_al2 = zeros(nsubs, nrois, nrois);
mst_pli_bet = zeros(nsubs, nrois, nrois);
mst_pli_gam = zeros(nsubs, nrois, nrois);

% --------------------construct mst-------------------------------------------------------
for sub = 1:nsubs
    mst_pli_del(sub,:,:) = kruskal_algorithm(squeeze(pli_delta_full_autofix(sub,:,:)));
    mst_pli_the(sub,:,:) = kruskal_algorithm(squeeze(pli_theta_full_autofix(sub,:,:)));
    mst_pli_al1(sub,:,:) = kruskal_algorithm(squeeze(pli_alpha1_full_autofix(sub,:,:)));
    mst_pli_al2(sub,:,:) = kruskal_algorithm(squeeze(pli_alpha2_full_autofix(sub,:,:)));
    mst_pli_bet(sub,:,:) = kruskal_algorithm(squeeze(pli_beta_full_autofix(sub,:,:)));
    mst_pli_gam(sub,:,:) = kruskal_algorithm(squeeze(pli_gamma_full_autofix(sub,:,:)));
end

nlrs = 6;
id = nrois:nrois:nrois*nlrs;

% ------------------pre-allocate supra-adjacency matrices---------------------------------
supra_mst_full = zeros(nrois*nlrs, nrois*nlrs, nsubs);
rand_supra_mst = zeros(nrois*nlrs, nrois*nlrs, nsubs);

% ------------------construct Minimum Spanning Tree (MST) supra-adjacency matrix ----------------------------------------
for c = 1:nsubs
    fprintf(1, 'Now constructing supra_mst for sub %s!\n', num2str(c))
    
    supratmp = blkdiag(squeeze(mst_pli_del(c,:,:)), squeeze(mst_pli_the(c,:,:)), ...
        squeeze(mst_pli_al1(c,:,:)), squeeze(mst_pli_al2(c,:,:)), squeeze(mst_pli_bet(c,:,:)), squeeze(mst_pli_gam(c,:,:)));
    for i = 1:length(id)
        supratmp(id(i)*size(supratmp,1)+1:size(supratmp,1)+1:end) = 1;
        supratmp(id(i)+1:size(supratmp, 1)+1:1+size(supratmp, 1)*min(size(supratmp, 1)-id(i),size(supratmp, 2))) = 1;
    end
    supra_mst_full(:,:,c) = supratmp;
end

 

% --------------------save mst supra-adjacency matrix-------------------------------------------
save('/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/supra_mst_full.mat','supra_mst_full') %import this into python.

%save raw pli matrices %collection of all 6 freq bands BL + FU in one
save('/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/full_raw_all.mat','*full_raw') %collection of all 6 freq bands BL + FU in one 


%%%%%% NOW -> RUN PYTHON SCRIPT TO CALCULATE MULTILAYER EIGENVECTOR CENTRALITY IN FPN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%% IMPORT OUTPUT PYTHON CSV FILE %%%%%%%%%%%%%%%%%%%
% Use script below to calculate average FPN eigenvector centrality per person per timepoint

% load multilayer EC
multilayer_ec_raw = readtable('/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/whole_multilayer_EC_mean.csv', 'HeaderLines', 1);

% place in SUBSxVALUES matrix %Values= amount of regions of FPN
id_ec = 1:12:889; %FPN: 1:12:889 12 fpn regios * 74 = 888 + 1 = 889. +1 because python starts at 0.
multilayer_ec_mat = zeros(74, 12); %(subjects, fpn_regions)
for i = 1:74
    multilayer_ec_mat(i,:) = multilayer_ec_raw{(id_ec(i):id_ec(i+1)-1), 2}';
end


save('/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/multilayer_nodal_ec.mat','multilayer_ec_mat')



% calculate average EC of fpn network
multilayer_ec_fpn_average = squeeze(mean(multilayer_ec_mat, 2));


save('/mnt/resource/m.vanlingen/m2b/matlab_output_final_okt2021/multilayer_ec_fpn_average.mat','multilayer_ec_fpn_average')



