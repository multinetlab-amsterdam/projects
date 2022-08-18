% Script to compute the correlations between simulations with average SC matrix 
% and individual empirical data for AEC and AEC full.

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

%% Script to correlate model FC ran with average SC to all individual correlations
% Loading model FC with average SC as input
% Loading of all empirical FCs - start with AEC over whole time series  
load('/path/to/data.mat')

%% Correlate AECf

coupling = 0.1:0.012:0.3;

for sub=1:numel(subs_emp)
    for cp = 1:numel(coupling)
        [rho_AECf(sub,cp), p_AECf(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR(cp,:,:))),vectorise(AECf_emp(:,:,sub)), 'type', 'spearman');
    end
end

for sub=1:numel(subs_emp)
    maxcor_AECf(sub)=max(rho_AECf(sub,:));
end

for sub=1:numel(subs_emp)
    index_p = find(rho_AECf(sub,:)==maxcor_AECf(sub));
    p_maxcor_AECf(sub) = p_AECf(sub,index_p); 
end 

for sub=1:numel(subs_emp)
    index_AECf(sub)=find(maxcor_AECf(sub)==rho_AECf(sub,:));
    opt_koppeling_AECf(sub) = coupling(index_AECf(sub)); 
end

%% Correlate AECe
for sub=1:numel(subs_emp)
    for cp = 1:numel(coupling)
        [rho_AECe(sub,cp), p_AECe(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR(cp,:,:))),vectorise(squeeze(AECe_emp(sub,:,:))), 'type', 'spearman');
    end
end

for sub=1:numel(subs_emp)
    maxcor_AECe(sub)=max(rho_AECe(sub,:));
end

for sub=1:numel(subs_emp)
    index_p = find(rho_AECe(sub,:)==maxcor_AECe(sub));
    p_maxcor_AECe(sub) = p_AECe(sub,index_p); 
end 

figure; 
plot(maxcor_AECf)
hold on 
plot(maxcor_AECe)

for sub=1:numel(subs_emp)
    index_AECe(sub)=find(maxcor_AECe(sub)==rho_AECe(sub,:));
    opt_koppeling_AECe(sub) = coupling(index_AECe(sub)); 
end

%% Test maximal correlations of individual SC runs vs average SC runs - for AECf and AECe

maxcor_AECe_av_SC = maxcor_AECe;
maxcor_AECf_av_SC = maxcor_AECf; 

[p_AECf,h_AECf,stats_AECf]=signrank(maxcor_AECf_av_SC,maxcor_AECf_N40);
[p_AECe,h_AECe,stats_AECe]=signrank(maxcor_AECe_av_SC,maxcor_AECe_N40);

figure; 
boxplot(maxcor_AECf_av_SC)
ylim([0.03 0.35])

figure;
boxplot(maxcor_AECf_N40)
ylim([0.03 0.35])

figure; 
boxplot(maxcor_AECe_av_SC)
ylim([0.05 0.3])

figure;
boxplot(maxcor_AECe_N40)
ylim([0.05 0.3])


%% Test whether there is a difference between AEC-full and AEC-e with paired test

[p_paired, h_paired, stats_paired]=signrank(maxcor_AECf, maxcor_AECe); 


