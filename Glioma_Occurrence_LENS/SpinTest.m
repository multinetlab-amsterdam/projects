function [pval, nullrho]= SpinTest(cortical, permno, atlas)
% SpinTest.m
% Computes the spin test for cortical data given a designated # of 
% permutations/spins of the input surface data in FreeSurfer fsaverage5
% space. The process is described in Alexander-Bloch et al 
% (2018, https://github.com/spin-test)
%
% FORMAT SpinTest(cortical, permno, atlas)
% cortical     - the csv filename of cortical data to spin (format 210x2).
%                  A csv file is expected.
% permno       - the number of permutations to generate the null
%                  distribution.
% atlas        - 'BNA' for brainnetome246 atlas, AAL for automatic 
%                  anatomical labeling atlas. Default: 'BNA' 
% 
% EXAMPLE: SpinTest('activity.csv', 100)
% will spin activity and tumor probability data, to 100 times from
% uniformly distributed sampled angles. 
%
%%author__ = Bernardo Maciel
%%contact__ = b.maciel@amsterdamumc.nl
%%date__ = 2021/11/22 
%%status__ = Production
%
%%%%%%%%%%%%%%%%%%%%
% Review History   %
%%%%%%%%%%%%%%%%%%%%
%
%%% Reviewed by Name Date
%
%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%
% Important: FREESURFER_HOME environment variable must be set.
% Expects ./atlas/ folder with freesurfer annotation and label files for
% both hemispheres for freesurfer and AAL atlases.
% Expects input data in ./data/ folder.
%%% Other m-files required:
% BNAcortical2surface.m
% SpinPermuFS.m (orginal from https://github.com/spin-test; it was adapted for this version and it's thus included in this project)
% nearestneighbour.m (https://nl.mathworks.com/matlabcentral/fileexchange/12574-nearestneighbour-m)

if ~exist('./analyses/', 'dir')
    mkdir('./analyses/')
end

filext = split(cortical, '.');
filext = split(filext(end-1), '/');

wsname = string(strcat('./analyses/tmp_', filext(end), '/'));
permname1 = strcat(wsname, 'rotation1.mat');
permname2 = strcat(wsname, 'rotation2.mat');

% Project cortical data to surface using BNA parcellation
cortical2surface(cortical, 1, wsname, atlas)

% Read the projected cortical data
surfaces = load(strcat(wsname,'data.mat'));
% Read the projected excluded vertices
exclusion = load(strcat(wsname,'exclude.mat'));

surface1_lh = surfaces.surface1_lh;  
surface1_rh = surfaces.surface1_rh; 
surface2_lh = surfaces.surface2_lh; 
surface2_rh = surfaces.surface2_rh;

% Spin surface data for both maps.
% Note: this requires the modified version of SpinPermuFS, the original
% one gets and loads filenames instead of variables to spin
SpinPermuFS(surface1_lh, surface1_rh, permno, permname1)
SpinPermuFS(surface2_lh, surface2_rh, permno, permname2)

surface1_lh(exclusion.x_left) = NaN; surface1_rh(exclusion.x_right) = NaN; 
surface2_lh(exclusion.x_left) = NaN; surface2_rh(exclusion.x_right) = NaN; 
% Load permutations. The name of the variables on the SpinPermuFS are
% bigrotr for permutations of the right hemisphere and bigrotl for the
% left
perm1 = load(permname1); perm1 = cat(2, perm1.bigrotl, perm1.bigrotr)';
perm2 = load(permname2); perm2 = cat(2, perm2.bigrotl, perm2.bigrotr)';

% Merge left and right surfaces into one
surface1 = cat(1, surface1_lh, surface1_rh); 
surface2 = cat(1, surface2_lh, surface2_rh); 

% Calculate the real Pearson's correlation between two maps
% 'rows','complete' to exclude NaN's
realrho = corr(surface1, surface2, 'type', 'spearman', 'rows', 'complete');

% Calculate correlation of the permutations of each map with the
% non-permuted value of the other. Overlap between the two maps is
% defined as the average of the two Pearson's correlation coefficients
nullrho1 = nan(permno, 1); nullrho2 = nan(permno, 1);

for i=1:permno
    nullrho1(i) = corr(perm1(:, i), surface2, 'type', 'spearman',  'rows', 'complete');
    nullrho2(i) = corr(perm2(:, i), surface1, 'type', 'spearman',  'rows', 'complete');
end

nullrho = mean([nullrho1 nullrho2], 2);
% Test the observed rho against null described by SpinPermuFS    
% assuming sign is preserved, calculate the probability that the observed
% correlation coeffcient is above the null distribution
pval = length(find(abs(nullrho)>abs(realrho)))/permno;
end
