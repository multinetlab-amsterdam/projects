% Script to run all correlations between simulations and empirical data.
% This script includes:
%        o       rho- and p-values of all correlations and maximum correlations
%        o       Median correlations per FC metric
%        o       The optimum coupling value per subject (coupling value that belongs to the maximum correlation)

%%%

%%author__ = Shanna Kulik & Preejas
%%contact__ = l.douw@amsterdamumc.nl
%%date__ = 2020
%%status__ = finished

%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

% Reviewed by Eduarda Centeno 2022/6

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

% Created and ran with Matlab R2018b

%%% Other m-files required:
% vectorise.m

%%% Toolboxes 
% None

%% Load data
% Concatenated loading model data AEC/PLI/PLV; 
% empirical data AEC; and empirical data PLI/PLV/AEC over epochs of 6 seconds
% and deleted unnecessary variables

load('/data.mat')

%% Corr FC
coupling = 0.1:0.012:0.3;

% limiting to 2 examples
mean_AEC_JR_example = mean_AEC_JR31(1:2,:,:,:);
mean_AEC_ep_emp_example = mean_AEC_ep_emp_31(1:2,:,:);
mean_AEC_full_emp_example = mean_AEC_full_31(1:2,:,:);
mean_PLI_JR_example = mean_PLI_JR31(1:2,:,:,:);
mean_PLI_ep_emp_example = mean_PLI_emp_31(1:2,:,:);
mean_PLV_JR_example = mean_PLV_JR31(1:2,:,:,:);
mean_PLV_ep_emp_example = mean_PLV_emp_31(1:2,:,:);

% correlation between NMM and empirical data Pearson 
for sub = 1:2
    for cp=1:numel(coupling)
        [rho_AEC_ep_example(sub,cp), p_AEC_ep_example(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR_example(sub,cp,:,:))),vectorise(squeeze(mean_AEC_ep_emp_example(sub,:,:))), 'type', 'Spearman');
        [rho_AEC_full_example(sub,cp), p_AEC_full_example(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR_example(sub,cp,:,:))),vectorise(squeeze(mean_AEC_full_emp_example(sub,:,:))), 'type', 'Spearman');
        [rho_PLI_example(sub,cp), p_PLI_example(sub,cp)] = corr(vectorise(squeeze(mean_PLI_JR_example(sub,cp,:,:))),vectorise(squeeze(mean_PLI_ep_emp_example(sub,:,:))), 'type', 'Spearman');
        [rho_PLV_example(sub,cp), p_PLV_example(sub,cp)]= corr(vectorise(squeeze(mean_PLV_JR_example(sub,cp,:,:))),vectorise(squeeze(mean_PLV_ep_emp_example(sub,:,:))), 'type', 'Spearman');
    end
end

for sub=1:2
    maxcor_AEC_ep_example(sub)=max(rho_AEC_ep_example(sub,:));
    maxcor_AEC_full_example(sub)=max(rho_AEC_full_example(sub,:));
    maxcor_PLI_example(sub)=max(rho_PLI_example(sub,:));
    maxcor_PLV_example(sub)=max(rho_PLV_example(sub,:));
end

for sub=1:2
    index_p_AEC_ep_example = find(rho_AEC_ep_example(sub,:)==maxcor_AEC_ep_example(sub));
    p_maxcor_AEC_ep_example(sub) = p_AEC_ep_example(sub,index_p_AEC_ep_example);
    index_p_AEC_full_example = find(rho_AEC_full_example(sub,:)==maxcor_AEC_full_example(sub));
    p_maxcor_AEC_full_example(sub) = p_AEC_full_example(sub,index_p_AEC_full_example);
    index_p_PLI_example = find(rho_PLI_example(sub,:)==maxcor_PLI_example(sub));
    p_maxcor_PLI_example(sub) = p_PLI_example(sub,index_p_PLI_example);
    index_p_PLV_example = find(rho_PLV_example(sub,:)==maxcor_PLV_example(sub));
    p_maxcor_PLV_example(sub) = p_PLV_example(sub,index_p_PLV_example);
end 


maxcor_AEC_ep_tot = maxcor_AEC_ep_example;
maxcor_AEC_full_tot = maxcor_AEC_full_example;
maxcor_PLI_tot = maxcor_PLI_example ;
maxcor_PLV_tot = maxcor_PLV_example ;

p_maxcor_AEC_ep_tot = p_maxcor_AEC_ep_example;
p_maxcor_AEC_full_tot = p_maxcor_AEC_full_example;
p_maxcor_PLI_tot = p_maxcor_PLI_example;
p_maxcor_PLV_tot = p_maxcor_PLV_example;

median_cor_AEC_ep_tot = median(maxcor_AEC_ep_tot); 
median_cor_AEC_full_tot = median(maxcor_AEC_full_tot); 
median_cor_PLI_tot = median(maxcor_PLI_tot); 
median_cor_PLV_tot = median(maxcor_PLV_tot); 
mean_AEC_ep_tot = mean(maxcor_AEC_ep_tot);


%% Optimal coupling parameters per sub for AEC long
for sub=1:2
    index_AEC_ep_example(sub)=find(maxcor_AEC_ep_example(sub)==rho_AEC_ep_example(sub,:));
    opt_koppeling_AEC_ep_example(sub) = coupling(index_AEC_ep_example(sub)); 
    index_AEC_full_example(sub)=find(maxcor_AEC_full_example(sub)==rho_AEC_full_example(sub,:));
    opt_koppeling_AEC_full_example(sub) = coupling(index_AEC_full_example(sub)); 
    index_PLI_example(sub)=find(maxcor_PLI_example(sub)==rho_PLI_example(sub,:));
    opt_koppeling_PLI_example(sub) = coupling(index_PLI_example(sub));
    index_PLV_example(sub)=find(maxcor_PLV_example(sub)==rho_PLV_example(sub,:));
    opt_koppeling_PLV_example(sub) = coupling(index_PLV_example(sub));
end
