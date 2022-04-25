% Brain activity and tumor occurrence

%   This script encompassess the analyses performed in the manuscript 
%   Regional healthy brain activity, glioma occurrence and symptomatology
%   by Numan et al. 

%%%

%%author__ = Tianne Numman & Linda Douw
%%contact__ = l.douw@amsterdamumc.nl
%%date__ = 2022/04/25 % date script was created
%%status__ = status % status of the script: Concluded/Finished


%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%

%%% Reviewed by 


%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%

%%% Other m-files required:
% spintest.m (also in this repo, based on the idea and m-files by Alexander-Bloch et al, see https://github.com/spin-test/spin-test)

%%% Toolboxes 
% Brain Connectivity Toolbox, version 2017_01_15
% (https://sites.google.com/site/bctnet/)

%%% Other
% We also use python scripts to calculate broadband power, offset and slope

%% First load all data

load('/data/KNW/l.douw/numan_lens/git_revision/LENS_source_data_final.mat')

%% Create intrinsic regional values using the HCs matched to the AMS cohort
    
% select relevant values 
HC_FOOOF_slope_matched_to_AMS = HC_FOOOF_slope(HC_matched_indx_AMS,:);
HC_FOOOF_offset_matched_to_AMS = HC_FOOOF_offset(HC_matched_indx_AMS,:);
HC_broadband_power_matched_to_AMS = HC_broadband_power(HC_matched_indx_AMS,:);

% create average of the 45 matched HCs
HC_FOOOF_slope_matched_to_AMS_mean = mean(HC_FOOOF_slope_matched_to_AMS,1);
HC_FOOOF_offset_matched_to_AMS_mean = mean(HC_FOOOF_offset_matched_to_AMS,1);
HC_broadband_power_matched_to_AMS_mean = mean(HC_broadband_power_matched_to_AMS,1);

%% Calculate LENS per individual patient in the AMS cohort  

nr_rois = 210;
ams_lens_bbp = nan(1,83);
ams_lens_offset = nan(1,83);
ams_lens_slope = nan(1,83);
nr_tumor_rois_AMS = nan(1,83);
mean_overlap_tumor_AMS = nan(1,83);
for i = 1:length(fname)
    overlap_summation = 0;
    for k = 1:length(fname(i,1).tumor_rois)
        if fname(i,1).tumor_rois(k) > nr_rois
        else
            fname(i,1).slope83(k) = HC_FOOOF_slope_matched_to_AMS_mean(fname(i,1).tumor_rois(k))*fname(i,1).perc_overlap(k);
            fname(i,1).offset83(k) = HC_FOOOF_offset_matched_to_AMS_mean(fname(i,1).tumor_rois(k))*fname(i,1).perc_overlap(k);
            fname(i,1).bbp83(k) = HC_broadband_power_matched_to_AMS_mean(fname(i,1).tumor_rois(k))*fname(i,1).perc_overlap(k);
            overlap_summation = overlap_summation + fname(i,1).perc_overlap(k); 
        end
    end
    ams_lens_slope(i) = nansum(fname(i,1).slope83)/overlap_summation;
    ams_lens_offset(i) = nansum(fname(i,1).offset83)/overlap_summation;
    ams_lens_bbp(i) = nansum(fname(i,1).bbp83)/overlap_summation;
    nr_tumor_rois_AMS(i) =  length(fname(i,1).tumor_rois);
    mean_overlap_tumor_AMS(i) = mean(fname(i,1).perc_overlap);
end

tumor_occ_AMS_perc = tumor_occ_AMS/83*100; % calculate percentage

[corr_slope_AMS, p_slope_AMS] = corr(tumor_occ_AMS_perc(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_AMS, p_offset_AMS] = corr(tumor_occ_AMS_perc(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_AMS, p_bbp_AMS] = corr(tumor_occ_AMS_perc(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');


%% subgroup analyses for AMS cohort

% subgroup 1: IDH mutant with codeletion
ams_data_idhmut_codel_cortical = ams_data_idhmut_codel(1:210,:);
ams_data_idhmut_codel_cortical_per = ams_data_idhmut_codel(1:210,:)/21*100;

% test associations with bbp/slope/offset
[corr_slope_ams_idhmut_codel, p_slope_ams_idhmut_codel] = corr(ams_data_idhmut_codel_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_ams_idhmut_codel, p_offset_ams_idhmut_codel] = corr(ams_data_idhmut_codel_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_ams_idhmut_codel, p_bbp_ams_idhmut_codel] = corr(ams_data_idhmut_codel_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

% subgroup 2: IDH mutant without codel
ams_data_idhmut_noncodel_cortical = ams_data_idhmut_noncodel(1:210,:);

% test associations with bbp/slope/offset
[corr_slope_ams_idhmut_noncodel, p_slope_ams_idhmut_noncodel] = corr(ams_data_idhmut_noncodel_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_ams_idhmut_noncodel, p_offset_ams_idhmut_noncodel] = corr(ams_data_idhmut_noncodel_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_ams_idhmut_noncodel, p_bbp_ams_idhmut_noncodel] = corr(ams_data_idhmut_noncodel_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

% subgroup 3: IDH wt
ams_data_idhwt_cortical = ams_data_idhwt(1:210,:);

% test associations with bbp/slope/offset
[corr_slope_ams_idhwt, p_slope_ams_idhwt] = corr(ams_data_idhwt_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_ams_idhwt, p_offset_ams_idhwt] = corr(ams_data_idhwt_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_ams_idhwt, p_bbp_ams_idhwt] = corr(ams_data_idhwt_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

%% ADDITIONAL ANALYSIS: does age impact any results for the significantly older subgroup 3?

HC_broadband_power_old = HC_broadband_power(age_hcs > 45,:);
HC_FOOOF_slope_old = HC_FOOOF_slope(age_hcs > 45,:);
HC_FOOOF_offset_old = HC_FOOOF_offset(age_hcs > 45,:);

[corr_slope_ams_idhwt_old, p_slope_ams_idhwt_old]  = corr(ams_data_idhwt_cortical(1:210), (mean(HC_FOOOF_slope_old,1))', 'type', 'spearman');
[corr_offset_ams_idhwt_old, p_offset_ams_idhwt_old] = corr(ams_data_idhwt_cortical(1:210), (mean(HC_FOOOF_offset_old,1))', 'type', 'spearman');
[corr_bbp_ams_idhwt_old, p_bbp_ams_idhwt_old]  = corr(ams_data_idhwt_cortical(1:210), (mean(HC_broadband_power_old,1))', 'type', 'spearman');

%% MGH cohort

boston_data_cortical = boston_data(1:210,:);
boston_data_cortical_per = boston_data_cortical/121*100; % calculate percentage

%% test associations of tumor occurrence with bbp/slope/offset

[corr_slope_boston, p_slope_boston] = corr(boston_data_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_boston, p_offset_boston] = corr(boston_data_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_boston, p_bbp_boston] = corr(boston_data_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

%% split by tumor subgroups

% subgroup 3: IDH wt
boston_data_idhwt_cortical = boston_data_idhwt(1:210,:);
boston_data_idhwt_cortical_per = boston_data_idhwt(1:210,:)/91*100;

% test associations with bbp/slope/offset
[corr_slope_boston_idhwt, p_slope_boston_idhwt] = corr(boston_data_idhwt_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_boston_idhwt, p_offset_boston_idhwt] = corr(boston_data_idhwt_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_boston_idhwt, p_bbp_boston_idhwt] = corr(boston_data_idhwt_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

% check with old HC cohort
[corr_slope_boston_idhwt_old, p_slope_boston_idhwt_old]  = corr(boston_data_idhwt_cortical(1:210), (mean(HC_FOOOF_slope_old,1))', 'type', 'spearman');
[corr_offset_boston_idhwt_old, p_offset_boston_idhwt_old] = corr(boston_data_idhwt_cortical(1:210), (mean(HC_FOOOF_offset_old,1))', 'type', 'spearman');
[corr_bbp_boston_idhwt_old, p_bbp_boston_idhwt_old]  = corr(boston_data_idhwt_cortical(1:210), (mean(HC_broadband_power_old,1))', 'type', 'spearman');

% subgroup 2: IDH mut noncodel
boston_data_idhmut_cortical = boston_data_idhmut(1:210,:);

% test associations with bbp/slope/offset
[corr_slope_boston_idhmut, p_slope_boston_idhmut] = corr(boston_data_idhmut_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_boston_idhmut, p_offset_boston_idhmut] = corr(boston_data_idhmut_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_boston_idhmut, p_bbp_boston_idhmut] = corr(boston_data_idhmut_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

%% calculate LENS for MGH data

nsubs_boston_idhwt = 91;
nr_rois = 210;
bos_tumor_load_bna_perc_idhwt = nan(nsubs_boston_idhwt,nr_rois);
for i = 1:nsubs_boston_idhwt
    for j = 1:nr_rois
        temp = bos_tumor_load_bna_idhwt(i,j)./bos_vols(j,1);
        bos_tumor_load_bna_perc_idhwt(i,j) = temp;
    end
end

bos_lens_bbp_idhwt = nan(nsubs_boston_idhwt,1);
bos_lens_offset_idhwt = nan(nsubs_boston_idhwt,1);
bos_lens_slope_idhwt = nan(nsubs_boston_idhwt,1);
for i = 1:nsubs_boston_idhwt
    bos_overlap_summation = 0;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhwt(i,:) .* HC_broadband_power_matched_to_AMS_mean(1,:)));
    bos_overlap_summation = bos_overlap_summation + nansum(bos_tumor_load_bna_perc_idhwt(i,:));
    bos_lens_bbp_idhwt(i,:) = temp/bos_overlap_summation;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhwt(i,:) .* HC_FOOOF_offset_matched_to_AMS_mean(1,:)));
    bos_lens_offset_idhwt(i,:) = temp/bos_overlap_summation;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhwt(i,:) .* HC_FOOOF_slope_matched_to_AMS_mean(1,:)));
    bos_lens_slope_idhwt(i,:) = temp/bos_overlap_summation;
end

%% calculate LENS for MGH data subgroup 2

nsubs_boston_idhmut = 18;
nr_rois = 210;
bos_tumor_load_bna_perc_idhmut = nan(nsubs_boston_idhmut,nr_rois);
for i = 1:nsubs_boston_idhmut
    for j = 1:nr_rois
        temp = bos_tumor_load_bna_idhmut(i,j)./bos_vols(j,1);
        bos_tumor_load_bna_perc_idhmut(i,j) = temp;
    end
end

bos_lens_bbp_idhmut = nan(nsubs_boston_idhmut,1);
bos_lens_offset_idhmut = nan(nsubs_boston_idhmut,1);
bos_lens_slope_idhmut = nan(nsubs_boston_idhmut,1);
for i = 1:nsubs_boston_idhmut
    bos_overlap_summation = 0;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhmut(i,:) .* HC_broadband_power_matched_to_AMS_mean(1,:)));
    bos_overlap_summation = bos_overlap_summation + nansum(bos_tumor_load_bna_perc_idhmut(i,:));
    bos_lens_bbp_idhmut(i,:) = temp/bos_overlap_summation;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhmut(i,:) .* HC_FOOOF_offset_matched_to_AMS_mean(1,:)));
    bos_lens_offset_idhmut(i,:) = temp/bos_overlap_summation;
    temp = nansum(nansum(bos_tumor_load_bna_perc_idhmut(i,:) .* HC_FOOOF_slope_matched_to_AMS_mean(1,:)));
    bos_lens_slope_idhmut(i,:) = temp/bos_overlap_summation;
end

%% TCGA-GBM cohort

tcga_gbm_data_cortical = tcga_data(1:210,:);
tcga_gbm_data_cortical_per = tcga_data/102*100; % calculate percentage

[corr_slope_tcga, p_slope_tcga] = corr(tcga_gbm_data_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_tcga, p_offset_tcga] = corr(tcga_gbm_data_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_tcga, p_bbp_tcga] = corr(tcga_gbm_data_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

%% TCGA-LGG cohort

tcga_data_idhmut_codel_cortical = tcga_data_idhmut_codel(1:210,:);
tcga_data_idhmut_noncodel_cortical = tcga_data_idhmut_noncodel(1:210,:);

% test associations with bbp/slope/offset
[corr_slope_tcga_idhmut_codel, p_slope_tcga_idhmut_codel] = corr(tcga_data_idhmut_codel_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman')
[corr_offset_tcga_idhmut_codel, p_offset_tcga_idhmut_codel] = corr(tcga_data_idhmut_codel_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman')
[corr_bbp_tcga_idhmut_codel, p_bbp_tcga_idhmut_codel] = corr(tcga_data_idhmut_codel_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman')


%% Do main analysis on the high SNR regions
% first find them

HC_bbp_q3 = quantile(HC_broadband_power_matched_to_AMS_mean,0.75);
HC_broadband_power_matched_mean_low = HC_broadband_power_matched_to_AMS_mean;
HC_broadband_power_matched_mean_low(HC_broadband_power_matched_to_AMS_mean > HC_bbp_q3) = NaN;
select_low_regions_bbp = find(HC_broadband_power_matched_to_AMS_mean <= HC_bbp_q3);

HC_slope_q3 = quantile(HC_FOOOF_slope_matched_to_AMS_mean,0.75);
HC_FOOOF_slope_matched_mean_low = HC_FOOOF_slope_matched_to_AMS_mean;
HC_FOOOF_slope_matched_mean_low(HC_FOOOF_slope_matched_to_AMS_mean > HC_slope_q3) = NaN;
select_low_regions_slope = find(HC_FOOOF_slope_matched_to_AMS_mean <= HC_slope_q3);

HC_offset_q3 = quantile(HC_FOOOF_offset_matched_to_AMS_mean,0.75);
HC_FOOOF_offset_matched_mean_low = HC_FOOOF_offset_matched_to_AMS_mean;
HC_FOOOF_offset_matched_mean_low(HC_FOOOF_offset_matched_to_AMS_mean > HC_offset_q3) = NaN;
select_low_regions_offset = find(HC_FOOOF_offset_matched_to_AMS_mean <= HC_offset_q3);

%%
% only q3
[corr_slope_tcga_idhmut_codel_q3, p_slope_tcga_idhmut_codel_q3] = corr(tcga_data_idhmut_codel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman')
[corr_offset_tcga_idhmut_codel_q3, p_offset_tcga_idhmut_codel_q3] = corr(tcga_data_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman')
[corr_bbp_tcga_idhmut_codel_q3, p_bbp_tcga_idhmut_codel_q3] = corr(tcga_data_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman')

% test associations with bbp/slope/offset
[corr_slope_tcga_idhmut_noncodel, p_slope_tcga_idhmut_noncodel] = corr(tcga_data_idhmut_noncodel_cortical(1:210), HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman')
[corr_offset_tcga_idhmut_noncodel, p_offset_tcga_idhmut_noncodel] = corr(tcga_data_idhmut_noncodel_cortical(1:210), HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman')
[corr_bbp_tcga_idhmut_noncodel, p_bbp_tcga_idhmut_noncodel] = corr(tcga_data_idhmut_noncodel_cortical(1:210), HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman')

% only q3
[corr_slope_tcga_idhmut_noncodel_q3, p_slope_tcga_idhmut_noncodel_q3] = corr(tcga_data_idhmut_noncodel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman')
[corr_offset_tcga_idhmut_noncodel_q3, p_offset_tcga_idhmut_noncodel_q3] = corr(tcga_data_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman')
[corr_bbp_tcga_idhmut_noncodel_q3, p_bbp_tcga_idhmut_noncodel_q3] = corr(tcga_data_idhmut_noncodel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman')


%% reproduce with old HCs

[corr_slope_tcga_old, p_slope_tcga_old] = corr(tcga_gbm_data_cortical(1:210), (mean(HC_FOOOF_slope_old,1))', 'type', 'spearman');
[corr_offset_tcga_old, p_offset_tcga_old] = corr(tcga_gbm_data_cortical(1:210), (mean(HC_FOOOF_offset_old,1))', 'type', 'spearman');
[corr_bbp_tcga_old, p_bbp_tcga_old] = corr(tcga_gbm_data_cortical(1:210), (mean(HC_broadband_power_old,1))', 'type', 'spearman');


%% Correlate tumor occurrence and brain activity measures with non-medial regions

% full cohorts
[corr_slope_ams_q3, p_slope_ams_q3] = corr(tumor_occ_AMS_perc(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_ams_q3, p_offset_ams_q3] = corr(tumor_occ_AMS_perc(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_ams_q3, p_bbp_ams_q3] = corr(tumor_occ_AMS_perc(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

[corr_slope_boston_q3, p_slope_boston_q3] = corr(boston_data_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_boston_q3, p_offset_boston_q3] = corr(boston_data_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_boston_q3, p_bbp_boston_q3] = corr(boston_data_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

[corr_slope_tcga_q3, p_slope_tcga_q3] = corr(tcga_gbm_data_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_tcga_q3, p_offset_tcga_q3] = corr(tcga_gbm_data_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_tcga_q3, p_bbp_tcga_q3] = corr(tcga_gbm_data_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

% subgroup 1
[corr_slope_ams_idhmut_codel_q3, p_slope_ams_idhmut_codel_q3] = corr(ams_data_idhmut_codel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_ams_idhmut_codel_q3, p_offset_ams_idhmut_codel_q3] = corr(ams_data_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_ams_idhmut_codel_q3, p_bbp_ams_idhmut_codel_q3] = corr(ams_data_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

% subgroup 2
[corr_slope_ams_idhmut_noncodel_q3, p_slope_ams_idhmut_noncodel_q3] = corr(ams_data_idhmut_noncodel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_ams_idhmut_noncodel_q3, p_offset_ams_idhmut_noncodel_q3] = corr(ams_data_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_ams_idhmut_noncodel_q3, p_bbp_ams_idhmut_noncodel_q3] = corr(ams_data_idhmut_noncodel_cortical(select_low_regions_bbp),(HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

[corr_slope_boston_idhmut_noncodel_q3, p_slope_boston_idhmut_noncodel_q3] = corr(boston_data_idhmut_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_boston_idhmut_noncodel_q3, p_offset_boston_idhmut_noncodel_q3] = corr(boston_data_idhmut_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_boston_idhmut_noncodel_q3, p_bbp_boston_idhmut_noncodel_q3] = corr(boston_data_idhmut_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

% subgroup 3
[corr_slope_ams_idhwt_q3, p_slope_ams_idhwt_q3] = corr(ams_data_idhwt_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_ams_idhwt_q3, p_offset_ams_idhwt_q3] = corr(ams_data_idhwt_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_ams_idhwt_q3, p_bbp_ams_idhwt_q3] = corr(ams_data_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

[corr_slope_boston_idhwt_q3, p_slope_boston_idhwt_q3] = corr(boston_data_idhwt_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_offset_boston_idhwt_q3, p_offset_boston_idhwt_q3] = corr(boston_data_idhwt_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_bbp_boston_idhwt_q3, p_bbp_boston_idhwt_q3] = corr(boston_data_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');


%% Investigate correlations between brain activity measures

exclude_high_regions_all = unique(cat(1,find(isnan(HC_broadband_power_matched_mean_low)),find(isnan(HC_FOOOF_offset_matched_mean_low)),find(isnan(HC_FOOOF_slope_matched_mean_low))));
select_low_regions_all = linspace(1,210,210);
select_low_regions_all(exclude_high_regions_all) = [];

[corr_bbp_slope_low,p_bbp_slope_low] = corr((HC_broadband_power_matched_to_AMS_mean(select_low_regions_all))',(HC_FOOOF_slope_matched_to_AMS_mean(select_low_regions_all))','type','spearman');
[corr_bbp_offset_low,p_bbp_offset_low] = corr((HC_broadband_power_matched_to_AMS_mean(select_low_regions_all))',(HC_FOOOF_offset_matched_to_AMS_mean(select_low_regions_all))','type','spearman');
[corr_offset_slope_low,p_offset_slope_low] = corr((HC_FOOOF_slope_matched_to_AMS_mean(select_low_regions_all))',(HC_FOOOF_offset_matched_to_AMS_mean(select_low_regions_all))','type','spearman');

%% Do correlations based on subgroups, not cohorts

all_idhwt_cortical = sum(cat(2,boston_data_idhwt_cortical,ams_data_idhwt_cortical,tcga_gbm_data_cortical),2);
all_idhmut_noncodel_cortical = sum(cat(2,ams_data_idhmut_noncodel_cortical,boston_data_idhmut_cortical,tcga_data_idhmut_noncodel_cortical),2);
all_idhmut_codel_cortical = sum(cat(2,ams_data_idhmut_codel_cortical,tcga_data_idhmut_codel_cortical),2);

[corr_slope_all_idhwt_q3, p_slope_all_idhwt_q3] = corr(all_idhwt_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_slope_all_idhmut_noncodel_q3, p_slope_all_idhmut_noncodel_q3] = corr(all_idhmut_noncodel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');
[corr_slope_all_idhmut_codel_q3, p_slope_all_idhmut_codel_q3] = corr(all_idhmut_codel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))', 'type', 'spearman');

[corr_offset_all_idhwt_q3, p_offset_all_idhwt_q3] = corr(all_idhwt_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_offset_all_idhmut_noncodel_q3, p_offset_all_idhmut_noncodel_q3] = corr(all_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');
[corr_offset_all_idhmut_codel_q3, p_offset_all_idhmut_codel_q3] = corr(all_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))', 'type', 'spearman');

[corr_bbp_all_idhwt_q3, p_bbp_all_idhwt_q3] = corr(all_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');
[corr_bbp_all_idhmut_noncodel_q3, p_bbp_all_idhmut_noncodel_q3] = corr(all_idhmut_noncodel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');
[corr_bbp_all_idhmut_codel_q3, p_bbp_all_idhmut_codel_q3] = corr(all_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))', 'type', 'spearman');

%% ADDITIONAL ANALYSIS: with all regions, including medial 

[corr_slope_all_idhwt, p_slope_all_idhwt] = corr(all_idhwt_cortical, HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_slope_all_idhmut_noncodel, p_slope_all_idhmut_noncodel] = corr(all_idhmut_noncodel_cortical, HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');
[corr_slope_all_idhmut_codel, p_slope_all_idhmut_codel] = corr(all_idhmut_codel_cortical, HC_FOOOF_slope_matched_to_AMS_mean', 'type', 'spearman');

[corr_offset_all_idhwt, p_offset_all_idhwt] = corr(all_idhwt_cortical, HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_all_idhmut_noncodel, p_offset_all_idhmut_noncodel] = corr(all_idhmut_noncodel_cortical, HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');
[corr_offset_all_idhmut_codel, p_offset_all_idhmut_codel] = corr(all_idhmut_codel_cortical, HC_FOOOF_offset_matched_to_AMS_mean', 'type', 'spearman');

[corr_bbp_all_idhwt, p_bbp_all_idhwt] = corr(all_idhwt_cortical, HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_all_idhmut_noncodel, p_bbp_all_idhmut_noncodel] = corr(all_idhmut_noncodel_cortical, HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');
[corr_bbp_all_idhmut_codel, p_bbp_all_idhmut_codel] = corr(all_idhmut_codel_cortical, HC_broadband_power_matched_to_AMS_mean', 'type', 'spearman');

%% Redo with old HC cohort in IDHwt

[corr_slope_all_idhwt_old, p_slope_all_idhwt_old] = corr(all_idhwt_cortical, (mean(HC_FOOOF_slope_old,1))', 'type', 'spearman');
[corr_offset_all_idhwt_old, p_offset_all_idhwt_old] = corr(all_idhwt_cortical, (mean(HC_FOOOF_offset_old,1))', 'type', 'spearman');
[corr_bbp_all_idhwt_old, p_bbp_all_idhwt_old] = corr(all_idhwt_cortical, (mean(HC_broadband_power_old,1))', 'type', 'spearman');

%% Create LENS for TCGA data

nsubs_tcga_gbm = 102;
nsubs_tcga_lgg = 107;

tcga_gbm_lens_bbp = nan(nsubs_tcga_gbm,1);
tcga_gbm_lens_offset = nan(nsubs_tcga_gbm,1);
tcga_gbm_lens_slope = nan(nsubs_tcga_gbm,1);
for i = 1:nsubs_tcga_gbm
    tcga_gbm_overlap_summation = 0;
    temp = nansum(nansum(tcga_gbm_tumor_load_bna_perc(i,:) .* HC_broadband_power_matched_to_AMS_mean(1,:)));
    tcga_gbm_overlap_summation = tcga_gbm_overlap_summation + nansum(tcga_gbm_tumor_load_bna_perc(i,:));
    tcga_gbm_lens_bbp(i,:) = temp/tcga_gbm_overlap_summation;
    temp = nansum(nansum(tcga_gbm_tumor_load_bna_perc(i,:) .* HC_FOOOF_offset_matched_to_AMS_mean(1,:)));
    tcga_gbm_lens_offset(i,:) = temp/tcga_gbm_overlap_summation;
    temp = nansum(nansum(tcga_gbm_tumor_load_bna_perc(i,:) .* HC_FOOOF_slope_matched_to_AMS_mean(1,:)));
    tcga_gbm_lens_slope(i,:) = temp/tcga_gbm_overlap_summation;
    
end

tcga_lgg_lens_bbp = nan(nsubs_tcga_lgg,1);
tcga_lgg_lens_offset = nan(nsubs_tcga_lgg,1);
tcga_lgg_lens_slope = nan(nsubs_tcga_lgg,1);
for i = 1:nsubs_tcga_lgg
    tcga_lgg_overlap_summation = 0;
    temp = nansum(nansum(tcga_lgg_tumor_load_bna_perc(i,:) .* HC_broadband_power_matched_to_AMS_mean(1,:)));
    tcga_lgg_overlap_summation = tcga_lgg_overlap_summation + nansum(tcga_lgg_tumor_load_bna_perc(i,:));
    tcga_lgg_lens_bbp(i,:) = temp/tcga_lgg_overlap_summation;
    temp = nansum(nansum(tcga_lgg_tumor_load_bna_perc(i,:) .* HC_FOOOF_offset_matched_to_AMS_mean(1,:)));
    tcga_lgg_lens_offset(i,:) = temp/tcga_lgg_overlap_summation;
    temp = nansum(nansum(tcga_lgg_tumor_load_bna_perc(i,:) .* HC_FOOOF_slope_matched_to_AMS_mean(1,:)));
    tcga_lgg_lens_slope(i,:) = temp/tcga_lgg_overlap_summation;
    
end

%% Create Figure 2A

% first calculate percentages

all_idhmut_codel_cortical_per = (all_idhmut_codel_cortical./48)*100;
all_idhmut_noncodel_cortical_per = (all_idhmut_noncodel_cortical./105)*100;
all_idhwt_cortical_per = (all_idhwt_cortical./207)*100;

% 

yneg = zeros(210,1);
xneg_bbp = nan(210,1);
xneg_offset = nan(210,1);
xneg_slope = nan(210,1);
for i = 1:210
    xneg_bbp(i,:) = std(HC_broadband_power_matched_to_AMS(:,i));
    xneg_offset(i,:) = std(HC_FOOOF_offset_matched_to_AMS(:,i));
    xneg_slope(i,:) = std(HC_FOOOF_slope_matched_to_AMS(:,i));
end


figure()
subplot(1,3,1)
scatter(all_idhmut_codel_cortical_per, HC_broadband_power_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_codel_cortical_per,HC_broadband_power_matched_to_AMS_mean,xneg_bbp,xneg_bbp,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,2)
scatter(all_idhmut_codel_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_codel_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,xneg_offset,xneg_offset,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,3)
scatter(all_idhmut_codel_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_codel_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,xneg_slope,xneg_slope,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create Figure 2B

subplot(1,3,1)
scatter(all_idhmut_noncodel_cortical_per, HC_broadband_power_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_noncodel_cortical_per,HC_broadband_power_matched_to_AMS_mean,xneg_bbp,xneg_bbp,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,2)
scatter(all_idhmut_noncodel_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_noncodel_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,xneg_offset,xneg_offset,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,3)
scatter(all_idhmut_noncodel_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhmut_noncodel_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,xneg_slope,xneg_slope,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create Figure 2C

subplot(1,3,1)
scatter(all_idhwt_cortical_per, HC_broadband_power_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhwt_cortical_per,HC_broadband_power_matched_to_AMS_mean,xneg_bbp,xneg_bbp,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,2)
scatter(all_idhwt_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhwt_cortical_per, HC_FOOOF_offset_matched_to_AMS_mean,xneg_offset,xneg_offset,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,3,3)
scatter(all_idhwt_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,'filled','k') 
hold all
errorbar(all_idhwt_cortical_per, HC_FOOOF_slope_matched_to_AMS_mean,xneg_slope,xneg_slope,yneg,yneg,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create Supplementary Figure 1

HC_broadband_power_matched_AMS_left = HC_broadband_power_matched_to_AMS(:,1:2:end);
HC_broadband_power_matched_AMS_right = HC_broadband_power_matched_to_AMS(:,2:2:end);

HC_FOOOF_offset_matched_AMS_left = HC_FOOOF_offset_matched_to_AMS(:,1:2:end);
HC_FOOOF_offset_matched_AMS_right = HC_FOOOF_offset_matched_to_AMS(:,2:2:end);

HC_FOOOF_slope_matched_AMS_left = HC_FOOOF_slope_matched_to_AMS(:,1:2:end);
HC_FOOOF_slope_matched_AMS_right = HC_FOOOF_slope_matched_to_AMS(:,2:2:end);

% determine errorbars per hemi

yneg_left = zeros(105,1);
xneg_bbp_left = nan(105,1);
xneg_offset_left = nan(105,1);
xneg_slope_left = nan(105,1);
yneg_right = zeros(105,1);
xneg_bbp_right = nan(105,1);
xneg_offset_right = nan(105,1);
xneg_slope_right = nan(105,1);
for i = 1:105
    xneg_bbp_left(i,:) = std(HC_broadband_power_matched_AMS_left(:,i));
    xneg_offset_left(i,:) = std(HC_FOOOF_offset_matched_AMS_left(:,i));
    xneg_slope_left(i,:) = std(HC_FOOOF_slope_matched_AMS_left(:,i));
    xneg_bbp_right(i,:) = std(HC_broadband_power_matched_AMS_right(:,i));
    xneg_offset_right(i,:) = std(HC_FOOOF_offset_matched_AMS_right(:,i));
    xneg_slope_right(i,:) = std(HC_FOOOF_slope_matched_AMS_right(:,i));
end

%% without q3

lowhigh_region_indices_bbp = zeros(210,1);
lowhigh_region_indices_bbp(select_low_regions_bbp) = 1;
lowhigh_region_indices_bbp_left = lowhigh_region_indices_bbp(1:2:end,:);
lowhigh_region_indices_bbp_right = lowhigh_region_indices_bbp(2:2:end,:);

lowhigh_region_indices_offset = zeros(210,1);
lowhigh_region_indices_offset(select_low_regions_offset) = 1;
lowhigh_region_indices_offset_left = lowhigh_region_indices_offset(1:2:end,:);
lowhigh_region_indices_offset_right = lowhigh_region_indices_offset(2:2:end,:);

lowhigh_region_indices_slope = zeros(210,1);
lowhigh_region_indices_slope(select_low_regions_slope) = 1;
lowhigh_region_indices_slope_left = lowhigh_region_indices_slope(1:2:end,:);
lowhigh_region_indices_slope_right = lowhigh_region_indices_slope(2:2:end,:);

%% Create Figure S1A

region_list = linspace(1,105,105);

subplot(1,2,1)
scatter(region_list, mean(HC_broadband_power_matched_AMS_left,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_bbp_left == 0), mean(HC_broadband_power_matched_AMS_left(:,lowhigh_region_indices_bbp_left == 0),1)','filled','r') 
errorbar(region_list,mean(HC_broadband_power_matched_AMS_left,1)',xneg_bbp_left,xneg_bbp_left,yneg_left,yneg_left,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,2,2)
scatter(region_list, mean(HC_broadband_power_matched_AMS_right,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_bbp_right == 0), mean(HC_broadband_power_matched_AMS_right(:,lowhigh_region_indices_bbp_right == 0),1)','filled','r') 
errorbar(region_list,mean(HC_broadband_power_matched_AMS_right,1)',xneg_bbp_right,xneg_bbp_right,yneg_right,yneg_right,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create Figure S1B

subplot(1,2,1)
scatter(region_list, mean(HC_FOOOF_offset_matched_AMS_left,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_offset_left == 0), mean(HC_FOOOF_offset_matched_AMS_left(:,lowhigh_region_indices_offset_left == 0),1)','filled','r') 
errorbar(region_list,mean(HC_FOOOF_offset_matched_AMS_left,1)',xneg_offset_left,xneg_offset_left,yneg_left,yneg_left,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,2,2)
scatter(region_list, mean(HC_FOOOF_offset_matched_AMS_right,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_offset_right == 0), mean(HC_FOOOF_offset_matched_AMS_right(:,lowhigh_region_indices_offset_right == 0),1)','filled','r') 
errorbar(region_list,mean(HC_FOOOF_offset_matched_AMS_right,1)',xneg_offset_right,xneg_offset_right,yneg_right,yneg_right,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create Figure S1C

subplot(1,2,1)
scatter(region_list, mean(HC_FOOOF_slope_matched_AMS_left,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_slope_left == 0), mean(HC_FOOOF_slope_matched_AMS_left(:,lowhigh_region_indices_slope_left == 0),1)','filled','r') 
errorbar(region_list,mean(HC_FOOOF_slope_matched_AMS_left,1)',xneg_slope_left,xneg_slope_left,yneg_left,yneg_left,'o','Color','#808080')
set(gca,'FontSize',20)

subplot(1,2,2)
scatter(region_list, mean(HC_FOOOF_slope_matched_AMS_right,1)','filled','k') 
hold all
scatter(region_list(lowhigh_region_indices_slope_right == 0), mean(HC_FOOOF_slope_matched_AMS_right(:,lowhigh_region_indices_slope_right == 0),1)','filled','r') 
errorbar(region_list,mean(HC_FOOOF_slope_matched_AMS_right,1)',xneg_slope_right,xneg_slope_right,yneg_right,yneg_right,'o','Color','#808080')
set(gca,'FontSize',20)

%% Create datasets to run spintest on 

%all regions
spin_slope_all_idhwt_q3 = cat(2,all_idhwt_cortical(select_low_regions_slope),(HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_slope_all_idhmut_codel_q3 = cat(2,all_idhmut_codel_cortical(select_low_regions_slope),(HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_slope_all_idhmut_noncodel_q3 = cat(2,all_idhmut_noncodel_cortical(select_low_regions_slope),(HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');

spin_offset_all_idhwt_q3 = cat(2,all_idhwt_cortical(select_low_regions_offset),(HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_offset_all_idhmut_codel_q3 = cat(2,all_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_offset_all_idhmut_noncodel_q3 = cat(2,all_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');

spin_bbp_all_idhwt_q3 = cat(2,all_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');
spin_bbp_all_idhmut_codel_q3 = cat(2,all_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');
spin_bbp_all_idhmut_noncodel_q3 = cat(2,all_idhmut_noncodel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

%AMS
spin_offset_ams_idhwt_q3 = cat(2,ams_data_idhwt_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_offset_ams_idhmut_codel_q3 = cat(2,ams_data_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_offset_ams_idhmut_noncodel_q3 = cat(2,ams_data_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');

spin_slope_ams_idhwt_q3 = cat(2,ams_data_idhwt_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_slope_ams_idhmut_codel_q3 = cat(2,ams_data_idhmut_codel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_slope_ams_idhmut_noncodel_q3 = cat(2,ams_data_idhmut_noncodel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');

spin_bbp_ams_idhwt_q3 = cat(2,ams_data_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');
spin_bbp_ams_idhmut_codel_q3 = cat(2,ams_data_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');
spin_bbp_ams_idhmut_noncodel_q3 = cat(2,ams_data_idhmut_noncodel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

%BOS
spin_offset_boston_idhwt_q3 = cat(2,boston_data_idhwt_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_offset_boston_idhmut_noncodel_q3 = cat(2,boston_data_idhmut_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');

spin_slope_boston_idhwt_q3 = cat(2,boston_data_idhwt_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_slope_boston_idhmut_noncodel_q3 = cat(2,boston_data_idhmut_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');

spin_bbp_boston_idhwt_q3 = cat(2,boston_data_idhwt_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');
spin_bbp_boston_idhmut_noncodel_q3 = cat(2,boston_data_idhmut_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

%TCGA-GBM
spin_slope_tcga_idhwt_q3 = cat(2,tcga_gbm_data_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_offset_tcga_idhwt_q3 = cat(2,tcga_gbm_data_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_bbp_tcga_idhwt_q3 = cat(2,tcga_gbm_data_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

%TCGA-LGG
spin_slope_tcga_idhmut_codel_q3 = cat(2,tcga_data_idhmut_codel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_offset_tcga_idhmut_codel_q3 = cat(2,tcga_data_idhmut_codel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_bbp_tcga_idhmut_codel_q3 = cat(2,tcga_data_idhmut_codel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

spin_slope_tcga_idhmut_noncodel_q3 = cat(2,tcga_data_idhmut_noncodel_cortical(select_low_regions_slope), (HC_FOOOF_slope_matched_mean_low(select_low_regions_slope))');
spin_offset_tcga_idhmut_noncodel_q3 = cat(2,tcga_data_idhmut_noncodel_cortical(select_low_regions_offset), (HC_FOOOF_offset_matched_mean_low(select_low_regions_offset))');
spin_bbp_tcga_idhmut_noncodel_q3 = cat(2,tcga_data_idhmut_noncodel_cortical(select_low_regions_bbp), (HC_broadband_power_matched_mean_low(select_low_regions_bbp))');

%% non q3

%all 
spin_slope_all_idhwt = cat(2,all_idhwt_cortical,HC_FOOOF_slope_matched_mean_83');
spin_slope_all_idhmut_codel = cat(2,all_idhmut_codel_cortical,HC_FOOOF_slope_matched_mean_83');
spin_slope_all_idhmut_noncodel = cat(2,all_idhmut_noncodel_cortical,HC_FOOOF_slope_matched_mean_83');

spin_offset_all_idhwt = cat(2,all_idhwt_cortical,(HC_FOOOF_offset_matched_mean_83)');
spin_offset_all_idhmut_codel = cat(2,all_idhmut_codel_cortical, (HC_FOOOF_offset_matched_mean_83)');
spin_offset_all_idhmut_noncodel = cat(2,all_idhmut_noncodel_cortical, (HC_FOOOF_offset_matched_mean_83)');

spin_bbp_all_idhwt = cat(2,all_idhwt_cortical, (HC_broadband_power_matched_mean_83)');
spin_bbp_all_idhmut_codel = cat(2,all_idhmut_codel_cortical, (HC_broadband_power_matched_mean_83)');
spin_bbp_all_idhmut_noncodel = cat(2,all_idhmut_noncodel_cortical, (HC_broadband_power_matched_mean_83)');

%AMS
spin_offset_ams_idhwt = cat(2,ams_data_idhwt_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_offset_ams_idhmut_codel = cat(2,ams_data_idhmut_codel_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_offset_ams_idhmut_noncodel = cat(2,ams_data_idhmut_noncodel_cortical, (HC_FOOOF_offset_matched_mean_AMS)');

spin_slope_ams_idhwt = cat(2,ams_data_idhwt_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_slope_ams_idhmut_codel = cat(2,ams_data_idhmut_codel_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_slope_ams_idhmut_noncodel = cat(2,ams_data_idhmut_noncodel_cortical, (HC_FOOOF_slope_matched_mean_AMS)');

spin_bbp_ams_idhwt = cat(2,ams_data_idhwt_cortical, (HC_broadband_power_matched_mean_AMS)');
spin_bbp_ams_idhmut_codel = cat(2,ams_data_idhmut_codel_cortical, (HC_broadband_power_matched_mean_AMS)');
spin_bbp_ams_idhmut_noncodel = cat(2,ams_data_idhmut_noncodel_cortical, (HC_broadband_power_matched_mean_AMS)');

%BOS
spin_offset_boston_idhwt = cat(2,boston_data_idhwt_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_offset_boston_idhmut_noncodel = cat(2,boston_data_idhmut_cortical, (HC_FOOOF_offset_matched_mean_AMS)');

spin_slope_boston_idhwt = cat(2,boston_data_idhwt_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_slope_boston_idhmut_noncodel = cat(2,boston_data_idhmut_cortical, (HC_FOOOF_slope_matched_mean_AMS)');

spin_bbp_boston_idhwt = cat(2,boston_data_idhwt_cortical, (HC_broadband_power_matched_mean_AMS)');
spin_bbp_boston_idhmut_noncodel = cat(2,boston_data_idhmut_cortical, (HC_broadband_power_matched_mean_AMS)');

%TCGA-GBM
spin_slope_tcga_idhwt = cat(2,tcga_gbm_data_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_offset_tcga_idhwt = cat(2,tcga_gbm_data_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_bbp_tcga_idhwt = cat(2,tcga_gbm_data_cortical, (HC_broadband_power_matched_mean_AMS)');

%TCGA-LGG
spin_slope_tcga_idhmut_codel = cat(2,tcga_data_idhmut_codel_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_offset_tcga_idhmut_codel = cat(2,tcga_data_idhmut_codel_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_bbp_tcga_idhmut_codel = cat(2,tcga_data_idhmut_codel_cortical, (HC_broadband_power_matched_mean_AMS)');

spin_slope_tcga_idhmut_noncodel = cat(2,tcga_data_idhmut_noncodel_cortical, (HC_FOOOF_slope_matched_mean_AMS)');
spin_offset_tcga_idhmut_noncodel = cat(2,tcga_data_idhmut_noncodel_cortical, (HC_FOOOF_offset_matched_mean_AMS)');
spin_bbp_tcga_idhmut_noncodel = cat(2,tcga_data_idhmut_noncodel_cortical, (HC_broadband_power_matched_mean_AMS)');

