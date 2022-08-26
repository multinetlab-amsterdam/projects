% Script to compute correlations with higher coupling from higher range for 
% all FC metrics

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
% vectorise.m

%%% Toolboxes 
% None

%% Load data DTI and MEG data

load('/path/to/data.mat')

for i=1:numel(subs_AECf)
    index(i) = find(strcmp({ASCIIS.name}, meg_AECf{i})==1);
end
AEC_emp = AEC(:,:,index); 

coupling = 0.2:0.012:0.4;

for sub = 1:numel(subs_AECf)
    for cp=1:numel(coupling)
        [rho_AECf(sub,cp), p_AECf(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR(sub,cp,:,:))),vectorise(AEC_emp(:,:,sub)), 'type', 'spearman');
    end
end

for sub=1:numel(subs_AECf)
    maxcor_AECf(sub)=max(rho_AECf(sub,:));
end

for sub=1:numel(subs_AECf)
    index_p = find(rho_AECf(sub,:)==maxcor_AECf(sub));
    p_maxcor_AECf(sub) = p_AECf(sub,index_p); 
end 

% find out whether the cases show a higher correlation between simulated 
% and empirical FC with a coupling value from the higher coupling range

index=find(maxcor_AECf(i)==rho_AECf(i,:));
opt_koppeling = coupling(index);