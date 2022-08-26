% Script to test between the maximum correlations of all FC metrics

%%%

%%author__ = Shanna Kulik & Prejaas Tewarie
%%contact__ = l.douw@amsterdamumc.nl
%%date__ = 2020
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2022/7

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% None

%%% Toolboxes 
% None

%% Script to test medians of all measures 

load('/path/to/data.mat')

%% Medians of each measure

med_AECf = median(maxcor_AECf_N40);
med_AECe = median(maxcor_AECe_N40);
med_PLI = median(maxcor_PLI_N40);
med_PLV = median(maxcor_PLV_N40); 

%% test medians of all measures or maximal correlations? 

[p_AECf_AECe,h_AECf_AECe,stats_AECf_AECe]=signrank(maxcor_AECf_N40,maxcor_AECe_N40);
[p_AECf_PLI,h_AECf_PLI,stats_AECf_PLI]=signrank(maxcor_AECf_N40,maxcor_PLI_N40);
[p_AECf_PLV,h_AECf_PLV,stats_AECf_PLV]=signrank(maxcor_AECf_N40,maxcor_PLV_N40);
[p_AECe_PLI,h_AECe_PLI,stats_AECe_PLI]=signrank(maxcor_AECe_N40,maxcor_PLI_N40);
[p_AECe_PLV,h_AECe_PLV,stats_AECe_PLV]=signrank(maxcor_AECe_N40,maxcor_PLV_N40);
[p_PLI_PLV,h_PLI_PLV,stats_PLI_PLV]=signrank(maxcor_PLI_N40,maxcor_PLV_N40);