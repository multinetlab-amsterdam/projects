% Script for shuffling data (AEC)

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
addpath('/path/Cross_sectional')
addpath('/path/Cross_disease')

% load('/path/Cross_sectional/AEC_theta_NO.mat')
load('/path/Cross_disease/AEC_HC_theta.mat')

AEC_HC_theta = AEC_HC_theta_tot; 
AEC_HC_theta(11:12,:,:,:,:)=AEC_theta_HC; 


% delete = [33 36 39 41 43 44 46 48 50 52 54 56 59 60 67 71 73 74 76 79]; 
% AEC_theta_NO(delete,:,:,:,:) = []; 
% AEC_theta_NO = (AEC_theta_NO+1)/2;

AEC_HC_theta = (AEC_theta_HC+1)/2;

for i=1:65
    for j=1:13
       for k=1:4
           AEC_HC_theta_zero(i,j,k,:,:)=squeeze(AEC_HC_theta(i,j,k,:,:)).*~eye(size(squeeze(AEC_HC_theta(i,j,k,:,:))));
       end 
    end
end
save('AEC_HC_theta_zero3', 'AEC_HC_theta_zero')

addpath('/mnt/anw-gold/KNW/s.kulik/MEG/HC')

% For NO: 
subj_NO = 1:71;
epochs_NO = 1:15;
% For HC:
subj_HC = 1:65;
epochs_HC = 1:13; 
sub_epochs = 1:4;

%% normalized clustering compare each epoch with the clustering of a 100 random networks with same edge weights.
nr_shuffle = 100; 
shuf_rand_matrix_cl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), nr_shuffle);
% shuf_rand_matrix_pl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), nr_shuffle);
cl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), 78);
% pl_NO = zeros(numel(subj_NO), numel(epochs_NO), numel(sub_epochs), 78); 
shuf_rand_matrix_cl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), nr_shuffle);
% shuf_rand_matrix_pl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), nr_shuffle);
cl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), 78);
% pl_HC = zeros(numel(subj_HC), numel(epochs_HC), numel(sub_epochs), 78); 


%% Loop for NO
tic 
for sub = 1:numel(subj_NO) %subs
    sub
    for ep = 1:numel(epochs_NO) %epochs
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(AEC_theta_NO_zero(sub,ep,sub_ep,:,:));
            for j = 1:nr_shuffle                                           %random shuffling 100x shuffle elements, shuffle in triu and then mirror shuffled triu to complete matrix
                up_ind = find(triu(ones(size(matrix,1)),1)==1);
                up_length = length(up_ind);
                up_elements = matrix(up_ind);
                matrix(up_ind) = up_elements(randperm(up_length));
                shuf_matrix = triu(matrix)+triu(matrix,1)';
%                calculate clustering on shuffled matrix
                shuf_rand_matrix_cl_NO(sub,ep,sub_ep,j) = mean(clustering_coef_wu(shuf_matrix));
%                 inv_matrix = 1./shuf_matrix;
%                 inv_matrix(~isfinite(inv_matrix))=0;
%                 D = distance_wei(inv_matrix);                               %the part that takes the longest
%                 shuf_rand_matrix_pl_NO(sub,ep,sub_ep,j) = charpath(D); %calculate diameter 
            end
            cl_NO(sub,ep,sub_ep,:) = mean(clustering_coef_wu(squeeze(AEC_theta_NO_zero(sub,ep,sub_ep,:,:))));
%             inv_m = 1./squeeze(AEC_theta_NO(sub,ep,sub_ep,:,:));
%             inv_m(~isfinite(inv_m))=0;
%             D2 = distance_wei(inv_m);
%             pl_NO(sub,ep,sub_ep,:) = charpath(D2);                
        end
    end
end
% 
% save('/path/Clus_Pl_NO_AEC_theta.mat', 'shuf_rand_matrix_cl_NO', 'shuf_rand_matrix_pl_NO', 'cl_NO', 'pl_NO')
save('/path/Clus_Pl_NO_AEC_theta_new_zero.mat', 'shuf_rand_matrix_cl_NO', 'cl_NO')


%% Loop for HC
% tic 
for sub = 1:numel(subj_HC) %subs
    sub
    for ep = 1:numel(epochs_HC) %epochs
        for sub_ep = 1:numel(sub_epochs)
            matrix = squeeze(AEC_HC_theta_zero(sub,ep,sub_ep,:,:));
            for j = 1:nr_shuffle                                           %random shuffling 100x shuffle elements, shuffle in triu and then mirror shuffled triu to complete matrix
                up_ind = find(triu(ones(size(matrix,1)),1)==1);
                up_length = length(up_ind);
                up_elements = matrix(up_ind);
                matrix(up_ind) = up_elements(randperm(up_length));
                shuf_matrix = triu(matrix)+triu(matrix,1)';
%                 figure; imagesc(shuf_matrix);colorbar; caxis([0 0.9])
%                 calculate clustering on shuffled matrix
                 shuf_rand_matrix_cl_HC(sub,ep,sub_ep,j) = mean(clustering_coef_wu(shuf_matrix));
%                 inv_matrix = 1./shuf_matrix;
%                 inv_matrix(~isfinite(inv_matrix))=0;
%                 D = distance_wei(inv_matrix);                               %the part that takes the longest
%                 shuf_rand_matrix_pl_HC(sub,ep,sub_ep,j) = charpath(D); %calculate diameter 
            end
            cl_HC(sub,ep,sub_ep,:) = mean(clustering_coef_wu(squeeze(AEC_HC_theta_zero(sub,ep,sub_ep,:,:))));
%             inv_m = 1./squeeze(AEC_HC_theta(sub,ep,sub_ep,:,:));
%             inv_m(~isfinite(inv_m))=0;
%             D2 = distance_wei(inv_m);
%             pl_HC(sub,ep,sub_ep,:) = charpath(D2);                
        end
    end
end

% % 
% % cl_HC(66:145,:,:,:) = []; 
% % pl_HC(66:145,:,:,:) = []; 
% % shuf_rand_matrix_cl_HC(66:145,:,:,:) = [];
% % shuf_rand_matrix_pl_HC(66:145,:,:,:) = [];
% 
% % save('/path/Clus_Pl_HC_AEC_theta.mat', 'shuf_rand_matrix_cl_HC', 'shuf_rand_matrix_pl_HC', 'cl_HC', 'pl_HC')
save('/path/Clus_HC_AEC3_theta_new_zero2020.mat', 'shuf_rand_matrix_cl_HC', 'cl_HC')


%% Calculate normalized CL and PL for pt and HC
% 
% norm_cl_HC = mean(cl_HC,4)./mean(shuf_rand_matrix_cl_HC,4);
% % norm_pl_HC = mean(pl_HC,4)./mean(shuf_rand_matrix_pl_HC,4);
% norm_cl_pt = mean(cl_NO,4)./mean(shuf_rand_matrix_cl_NO,4);
% % norm_pl_pt = mean(pl_NO,4)./mean(shuf_rand_matrix_pl_NO,4);
% % 
% save('/path/norm_HC_theta2.mat', 'norm_cl_HC')
% save('/path/norm_NO_theta2.mat', 'norm_cl_pt')
