% -------------------------------------------------------------
% tumor_tract_density_map_indices
% Script to compute Lesion-Tract Density Map (L-TDM)
% and Lesion-Tract Density Index (L-TDI)
% using only tumor-intersecting streamlines.
%
% Input: fibers_vox (streamline voxel coords, cell array)
%        fibers_wt (weights >0 means fiber intersects tumor)
%        tumor mask (NIfTI, already loaded as tumor_mask)
% -------------------------------------------------------------
%%%

%%author__ = Mona Zimmermann
%%contact__ = m.l.m.zimmermann@amsterdamumc.nl
%%date__ = 2026/04/30
%%status__ = Finished 

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%
%%% spm

%% --- USER EDIT: paths ---

function tumor_tract_density_map_indices_main_for_sbatch(path, sub_id) 

    load('/path/to/dTOR_fibers_vox_half_mm.mat'); % Provide path and filename for voxel-resampled connectome file (this determines the resolution in which the streamline sampling will be conducted)
    
    
    fiber_wt_mat   = fullfile(path, sub_id, [sub_id, '_tumormask2MNI152NLin2009bAsym_res-05mm_fiber_wts.mat']); % contains fibers_wt
    
    tumor_mask_file = fullfile(path,sub_id, [sub_id,'_tumormask2MNI152NLin2009bAsym_res-05mm.nii.gz']); % Insert ROI file here 
    
    %unzip .nii.gz tumor mask file
    gunzip(tumor_mask_file);
    tumor_mask_file_unzipped = strrep(tumor_mask_file,'.gz','');
    
    load(fiber_wt_mat, 'fibers_wt', 'fibers_vox'); % make sure both variables exist in file
    
    v = spm_vol(tumor_mask_file_unzipped);          % header
    mask_img = spm_read_vols(v);     % actual 3D mask
    tumor_mask = mask_img > 0;       % binarize just to be sure (tumor mask is already binarized)
    
    
    % Preallocate lesion-specific map
    L_TDM = zeros(size(tumor_mask));
    
    disp('Building lesion-specific tract density map...');
    
    for i = 1:length(fibers_vox)
        if fibers_wt(i) > 0   % only include tumor-intersecting streamlines
            streamline = fibers_vox{i};
    
            % keep only in-bounds voxels
            valid_idx = all(streamline > 0, 2) & ...
                        streamline(:,1) <= size(L_TDM,1) & ...
                        streamline(:,2) <= size(L_TDM,2) & ...
                        streamline(:,3) <= size(L_TDM,3);
            streamline = streamline(valid_idx,:);
    
            if isempty(streamline), continue; end
    
            % linear indices for voxel accumulation
            lin_idx = sub2ind(size(L_TDM), ...
                              streamline(:,1), ...
                              streamline(:,2), ...
                              streamline(:,3));
            L_TDM(lin_idx) = L_TDM(lin_idx) + 1;
        end
    
        if rem(i,10000)==0
            disp(['Processed streamline: ', num2str(i)]);
        end
    end
    
    disp('Lesion-specific density map built.');
    
    % Compute L-TDI (average lesion-specific density across the brain)
    %lesion_vox = L_TDM(tumor_mask > 0);
    %lesion_vox = L_TDM(L_TDM > 0);
    L_TDI = mean(L_TDM(L_TDM > 0));
    
    disp(['L-TDI = ', num2str(L_TDI)]);
    
    % Save outputs
    outdir = path;
    %[~, name, ~] = fileparts(tumor_mask_file_unzipped);
    
    Vout = v; % reuse header from mask
    Vout.fname = fullfile(outdir, sub_id, [sub_id '_L_TDM_MNI152NLin2009bAsym_res-05mm.nii']);
    spm_write_vol(Vout, L_TDM);
    
    save(fullfile(outdir, sub_id, [sub_id '_L_TDI_MNI152NLin2009bAsym_res-05mm.mat']), 'L_TDI');

    fprintf('Yaaaaay subject is finished %s:\n', sub_id);
