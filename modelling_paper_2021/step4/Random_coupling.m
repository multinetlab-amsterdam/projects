% Script to run comparison between maximum correlations of a subject’s own 
% data with maximum correlations of simulated data with empirical data of 
% other subjects. Also includes the figure (forest plot) that represents this (figure 5 in the manuscript).

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

%% Load all model FC's of AEC

load('/path/to/data.mat')


%% Get maximal coupling values for AECf and AECe
coupling = 0.1:0.012:0.3; 

for sub = 1:40
    for cp=1:numel(coupling)
        [rho_AECf(sub,cp), p_AECf(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR_tot(sub,cp,:,:))),vectorise(squeeze(AECf_emp(:,:,sub))), 'type', 'Spearman');
        [rho_AECe(sub,cp), p_AECe(sub,cp)] = corr(vectorise(squeeze(mean_AEC_JR_tot(sub,cp,:,:))),vectorise(squeeze(mean_AECe_emp(sub,:,:))), 'type', 'Spearman');
    end
end

for sub=1:40
    maxcor_AECf_tot(sub)=max(rho_AECf(sub,:));
    maxcor_AECe_tot(sub)=max(rho_AECe(sub,:));
end

for sub=1:40
    index_AECf(sub)=find(maxcor_AECf_tot(sub)==rho_AECf(sub,:));
    opt_koppeling_AECf_tot(sub) = coupling(index_AECf(sub)); 
    index_AECe(sub)=find(maxcor_AECe_tot(sub)==rho_AECe(sub,:));
    opt_koppeling_AECe_tot(sub) = coupling(index_AECe(sub));
end

%% Match FC simulations randomly

list = 1:40; 

% % Get model FC matrices with optimal coupling value 
for i=1:40
    mean_AECf_JR_opt(i,:,:) = squeeze(mean_AEC_JR_tot(i,index_AECf(i),:,:)); 
    mean_AECe_JR_opt(i,:,:) = squeeze(mean_AEC_JR_tot(i,index_AECe(i),:,:));
end

%% Do correlations for same participant but also for all other participants

for sub=1:40
    sub
    for sub2=1:40
        [rho_AECf_rand(sub,sub2), p_AECf_rand(sub,sub2)] = corr(vectorise(squeeze(mean_AECf_JR_opt(sub,:,:))), vectorise(squeeze(AECf_emp(:,:,sub2))), 'type', 'Spearman');
        [rho_AECe_rand(sub,sub2), p_AECe_rand(sub,sub2)] = corr(vectorise(squeeze(mean_AECe_JR_opt(sub,:,:))), vectorise(squeeze(mean_AECe_emp(sub2,:,:))), 'type', 'Spearman');
    end
end

% The diagonal in these matrices is the correlation with the own empirical
% data

%% Compare per subject the random correlations with the individual correlations

% Rank correlations with random other persons 
for i=1:40
    rank_rho_AECf_rand(i,:) = sort(rho_AECf_rand(i,:));
    rank_rho_AECe_rand(i,:) = sort(rho_AECe_rand(i,:));
end

% Self correlations must be in the top 97.5% to be significant - out of 40
% values ​​is the top 97.5% from the 39th value

for sub=1:40
    rank_AECf(sub)=find(rank_rho_AECf_rand(sub,:) == maxcor_AECf_tot(sub)); 
    rank_AECe(sub)=find(rank_rho_AECe_rand(sub,:) == maxcor_AECe_tot(sub));
end 

test_AECf = find(rank_AECf>38); % 5 out of 40 
test_AECe = find(rank_AECe>38); % 5 out of 40 

%% Make figure - forest plot - AECf

% What do we need:
% max and min correlation per person: horizontal bar
% median correlation: square or something
% Own correlation: red asterisk

max_cor_sub_AECf = max(rho_AECf_rand, [], 2); 
min_cor_sub_AECf = min(rho_AECf_rand, [],2); 

for i=1:40
    dis_min_max_AECf(i) = max_cor_sub_AECf(i) - min_cor_sub_AECf(i);      % Distance between minimum and maximum correlation per person
    med_AECf(i) = median(rho_AECf_rand(i,:));                             % Median 
    dis_real_med(i) = rho_AECf_rand(i,i)-med_AECf(i);                     % Distance between own correlations and median
end

for i=1:40
    for j=1:40
        rho_AECf_rand_min_med(i,j) = rho_AECf_rand(i,j)-med_AECf(i);      % Distance for each value between correlations and median per person
    end
end

for i=1:40
    sort_rho_AECf_rand(i,:) = sort(rho_AECf_rand_min_med(i,:));           % Sort values ​​of all correlations minus the median
end 
rho_AECf_rand_95 = sort_rho_AECf_rand(:,2:39);                            % Subtract 1st and last value

max_cor_sub_AECf_95 = max(rho_AECf_rand_95, [],2); 
min_cor_sub_AECf_95 = min(rho_AECf_rand_95, [],2); 

[sort_dis_real_med, index_dis_real_med] = sort(dis_real_med);             % Sorting distance between own correlations and median - greatest distance first

diag_rho_AECf_rand = diag(rho_AECf_rand);
for i=1:40
    diag_rho_AECf_min_med(i) = diag_rho_AECf_rand(i)-med_AECf(i); 
end

sort_diag_rho_AECf_min_med = diag_rho_AECf_min_med(index_dis_real_med);
sort_min_cor_sub_AECf_95 = min_cor_sub_AECf_95(index_dis_real_med);
sort_max_cor_sub_AECf_95 = max_cor_sub_AECf_95(index_dis_real_med);
sort_med_rho_AECf_med= med_AECf(index_dis_real_med); 

figure;
for i=1:40   
    plot([sort_min_cor_sub_AECf_95(i) sort_max_cor_sub_AECf_95(i)],[i,i], 'Color', [0.25 0.25 0.25], 'LineWidth', 2)
    hold on
    plot(0,i,'sk', 'MarkerSize', 10,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0 0 0]) 
    plot(sort_diag_rho_AECf_min_med(i),i,'rd', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0])
    ylim([0 41])
    xlim([-0.16 0.2])
    xticks([-0.16 0 0.2])
    xticklabels({'low', 'moderate', 'high'})
    yticks([1 40])
    yticklabels({'min. specificity', 'max. specificity'})
end

width = 8;     % Width in inches
height = 7;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 13;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

%% Make figure - forest plot - AECe

max_cor_sub_AECe = max(rho_AECe_rand, [], 2); 
min_cor_sub_AECe = min(rho_AECe_rand, [],2); 

for i=1:40
    dis_min_max_AECe(i) = max_cor_sub_AECe(i) - min_cor_sub_AECe(i);       % Distance between minimum and maximum correlation per person
    med_AECe(i) = median(rho_AECe_rand(i,:));                              % Median 
    dis_real_med_AECe(i) = rho_AECe_rand(i,i)-med_AECe(i);                 % Distance between own correlations and median
end

for i=1:40
    for j=1:40
        rho_AECe_rand_min_med(i,j) = rho_AECe_rand(i,j)-med_AECe(i);       % Distance for each value between correlations and median per person
    end
end

for i=1:40
    sort_rho_AECe_rand(i,:) = sort(rho_AECe_rand_min_med(i,:));            % Sort values ​​of all correlations minus the median
end 
rho_AECe_rand_95 = sort_rho_AECe_rand(:,2:39);                             % Subtract 1st and last value

max_cor_sub_AECe_95 = max(rho_AECe_rand_95, [],2); 
min_cor_sub_AECe_95 = min(rho_AECe_rand_95, [],2); 

[sort_dis_real_med_AECe, index_dis_real_med_AECe] = sort(dis_real_med_AECe);   % Sort distance between own correlations and median - greatest distance first

diag_rho_AECe_rand = diag(rho_AECe_rand);
for i=1:40
    diag_rho_AECe_min_med(i) = diag_rho_AECe_rand(i)-med_AECe(i); 
end

sort_diag_rho_AECe_min_med = diag_rho_AECe_min_med(index_dis_real_med_AECe);
sort_min_cor_sub_AECe_95 = min_cor_sub_AECe_95(index_dis_real_med_AECe);
sort_max_cor_sub_AECe_95 = max_cor_sub_AECe_95(index_dis_real_med_AECe);
sort_med_rho_AECe_med= med_AECf(index_dis_real_med_AECe); 

figure;
for i=1:40
    plot([sort_min_cor_sub_AECe_95(i) sort_max_cor_sub_AECe_95(i)],[i,i], 'Color', [0.25 0.25 0.25], 'LineWidth', 2)
    hold on
    plot(0,i,'sk', 'MarkerSize', 10,'MarkerEdgeColor', 'black', 'MarkerFaceColor', [0 0 0])
    plot(sort_diag_rho_AECe_min_med(i),i,'rd', 'MarkerSize', 10, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 0])
    ylim([0 41])
    xticks([-0.15 0 0.15])
    xticklabels({'low', 'moderate', 'high'})
    yticks([1 40])
    yticklabels({'min. specificity', 'max. specificity'})

end


width = 8;     % Width in inches
height = 7;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 13;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);