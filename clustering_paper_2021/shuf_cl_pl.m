% Script for shuffling data (PLI)

%   Calculate the weighted path length and clustering of the average HC
%   matrix - reshuffle matrix to normalize

%%%

%%author__ = Shanna Kulik & Jolanda Derks
%%contact__ = s.kulik@amsterdamumc.nl
%%date__ = 2021/3/7
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2021/11/4

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% NA

%%% Toolboxes 
% Brain Connectivity Toolbox https://sites.google.com/site/bctnet/ 
% version 2017-15-01


%% Clean initialization
clear 
close all

%% Add local path for BCT
addpath('/path/bct')

%% Load data 
addpath('/path/')
% load('/path/PLIallnew.mat')

addpath('/path/HC/')
load('/path/HC/PLI_theta_HC_new_2020.mat')

% PLI_NO_alpha = PLI_alpha; 
% % PLI_NO_theta = PLI_theta;
% delete = [33 36 39 41 43 44 46 48 50 52 54 56 59 60 67 71 73 74 76 79]; 

% PLI_NO_theta(delete,:,:,:,:) = []; 
% PLI_NO_alpha(delete,:,:,:,:) = []; 

% PLI_HC_alpha = PLI_alpha;

PLI_HC_theta2 = PLI_HC_theta(11:12,:,:,:,:);

% For NO: 
% subj_NO = 1:71;
% epochs_NO = 1:15;
% For HC:
subj_HC = 1:65;
epochs_HC = 1:13; 
sub_epochs = 1:4;

%% normalized clustering compare each epoch with the clustering of a 100 random networks with same edge weights.
nr_shuffle = 100; 
% shuf_rand_matrix_cl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), nr_shuffle);
% shuf_rand_matrix_pl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), nr_shuffle);
% cl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), 78);
% pl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), 78); 
shuf_rand_matrix_cl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), nr_shuffle);
shuf_rand_matrix_pl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), nr_shuffle);
cl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), 78);
pl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), 78); 

%% Loop for NO
tic 
for sub = 1:numel(subj_NO) %subs
    sub
    for ep = 1:numel(epochs_NO) %epochs
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(PLI_NO_alpha(sub,ep,sub_ep,:,:));
            for j = 1:nr_shuffle                                           %random shuffling 100x shuffle elements, shuffle in triu and then mirror shuffled triu to complete matrix
                up_ind = find(triu(ones(size(matrix,1)),1)==1);
                up_length = length(up_ind);
                up_elements = matrix(up_ind);
                matrix(up_ind) = up_elements(randperm(up_length));
                shuf_matrix = triu(matrix)+triu(matrix,1)';
%                calculate clustering on shuffled matrix
                shuf_rand_matrix_cl_NO(sub,ep,sub_ep,j) = mean(clustering_coef_wu(shuf_matrix));
                inv_matrix = 1./shuf_matrix;
                inv_matrix(~isfinite(inv_matrix))=0;
                D = distance_wei(inv_matrix);                               %takes the longest in this loop
                shuf_rand_matrix_pl_NO(sub,ep,sub_ep,j) = charpath(D); %calculate diameter 
            end
            cl_NO(sub,ep,sub_ep,:) = mean(clustering_coef_wu(squeeze(PLI_NO_alpha(sub,ep,sub_ep,:,:))));
            inv_m = 1./squeeze(PLI_NO_alpha(sub,ep,sub_ep,:,:));
            inv_m(~isfinite(inv_m))=0;
            D2 = distance_wei(inv_m);
            pl_NO(sub,ep,sub_ep,:) = charpath(D2);                
        end
    end
end

% 
% save('/path/Clus_Pl_NO_alpha.mat', 'shuf_rand_matrix_cl_NO', 'shuf_rand_matrix_pl_NO', 'cl_NO', 'pl_NO')

%% Loop for HC
tic 
for sub = 1 %subs
    sub
    for ep = 1:numel(epochs_HC) %epochs
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(PLI_HC_theta(sub,ep,sub_ep,:,:));
            figure; imagesc(matrix);colorbar; caxis([0 0.9])
            for j = 1:10 %nr_shuffle                                           %random shuffling 100x shuffle elements, shuffle in triu and then mirror shuffled triu to complete matrix
                [R,eff]=randmio_und(matrix, 1000);
                up_ind = find(triu(ones(size(matrix,1)),1)==1);
                up_length = length(up_ind);
                up_elements = matrix(up_ind);
                matrix(up_ind) = up_elements(randperm(up_length));
                shuf_matrix = triu(matrix)+triu(matrix,1)';
%                 figure; imagesc(R);colorbar; caxis([0 0.9])
%                 calculate clustering on shuffled matrix
                shuf_rand_matrix_cl_HC(sub,ep,sub_ep,j) = mean(clustering_coef_wu(shuf_matrix));
                inv_matrix = 1./shuf_matrix;
                inv_matrix(~isfinite(inv_matrix))=0;
                D = distance_wei(inv_matrix);                               %takes the longest in this loop
                shuf_rand_matrix_pl_HC(sub,ep,sub_ep,j) = charpath(D); %calculate diameter 
            end
            cl_HC(sub,ep,sub_ep,:) = mean(clustering_coef_wu(squeeze(PLI_HC_theta2(sub,ep,sub_ep,:,:))));
            inv_m = 1./squeeze(PLI_HC_theta2(sub,ep,sub_ep,:,:));
            inv_m(~isfinite(inv_m))=0;
            D2 = distance_wei(inv_m);
            pl_HC(sub,ep,sub_ep,:) = charpath(D2);                
        end
    end
end

% 
% cl_HC(66:145,:,:,:) = []; 
% pl_HC(66:145,:,:,:) = []; 
% shuf_rand_matrix_cl_HC(66:145,:,:,:) = [];
% shuf_rand_matrix_pl_HC(66:145,:,:,:) = [];

save('/path/Clus_Pl_HC_theta_new_2020.mat', 'shuf_rand_matrix_cl_HC', 'shuf_rand_matrix_pl_HC', 'cl_HC', 'pl_HC')


%% Calculate normalized CL and PL for pt and HC
% 
% % norm_cl_HC = mean(cl_HC,4)./mean(shuf_rand_matrix_cl_HC,4);
% % norm_pl_HC = mean(pl_HC,4)./mean(shuf_rand_matrix_pl_HC,4);
% norm_cl_pt = mean(cl_NO,4)./mean(shuf_rand_matrix_cl_NO,4);
% norm_pl_pt = mean(pl_NO,4)./mean(shuf_rand_matrix_pl_NO,4);
% 
% % save('/path/norm_HC_theta.mat', 'norm_cl_HC', 'norm_pl_HC')
% save('/path/norm_NO_alpha.mat', 'norm_cl_pt','norm_pl_pt')
