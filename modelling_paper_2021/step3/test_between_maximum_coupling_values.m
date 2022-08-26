% Script to apply friedman test and check whether the maximum coupling 
% values per FC metric differ

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

%% Test whether there is a difference between coupling values belonging to the maximum correlations

load('/path/to/data.mat')

opt_koppeling_all(:,1)=opt_koppeling_AECf_tot; 
opt_koppeling_all(:,2)=opt_koppeling_AECe_tot; 
opt_koppeling_all(:,3)=opt_koppeling_PLI_tot; 
opt_koppeling_all(:,4)=opt_koppeling_PLV_tot; 


[p, tbl, stats]=friedman(opt_koppeling_all,4); 
