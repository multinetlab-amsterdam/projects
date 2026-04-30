%% --- Load the tractogram ---
load('/path/to/dTOR_fibers_vox_half_mm.mat', 'fibers_vox');

%% --- Load MNI template to get size and header info ---
mni_file = '/path/to/mni_icbm152_t1_tal_nlin_asym_09b_hires.nii'; % adjust path to your NIfTI
mni_info = niftiinfo(mni_file);
template_size = mni_info.ImageSize; % automatically get template size

%% --- Initialize Tract Density Map ---
TDM = zeros(template_size, 'uint32'); % uint32 to save memory

%% --- Populate TDM ---
num_fibers = numel(fibers_vox);
for f = 1:num_fibers
    voxels = fibers_vox{f}; % Nx3 array of x/y/z coordinates
    
    % Round coordinates to nearest voxel
    voxels = round(voxels);
    
    % Keep only valid indices within template
    valid_idx = all(voxels >= 1, 2) & ...
                voxels(:,1) <= template_size(1) & ...
                voxels(:,2) <= template_size(2) & ...
                voxels(:,3) <= template_size(3);
    voxels = voxels(valid_idx, :);
    
    % Convert 3D indices to linear indices
    lin_idx = sub2ind(template_size, voxels(:,1), voxels(:,2), voxels(:,3));
    
    % Increment TDM
    TDM(lin_idx) = TDM(lin_idx) + 1;
    
    % Optional progress display
    if mod(f, 100000) == 0
        fprintf('Processed %d/%d fibers...\n', f, num_fibers);
    end
end

%% --- Save TDM as a NIfTI ---
tdm_file = '/path/to/whole_brain_TDM_mni_icbm152_nlin_asym_09b_new_script_final.nii';
% Copy header info from MNI template
tdm_info = mni_info;
tdm_info.Datatype = 'uint32';
tdm_info.BitDepth = 32;
tdm_info.ImageSize = size(TDM);
tdm_info.PixelDimensions = [0.5 0.5 0.5];

% Kill intensity scaling
tdm_info.AdditiveOffset = 0;
tdm_info.MultiplicativeScaling = 1;
tdm_info.raw.scl_inter = 0;
tdm_info.raw.scl_slope = 1;


niftiwrite(TDM, tdm_file, tdm_info,'Compressed', true);
fprintf('Tract density map saved as NIfTI: %s\n', tdm_file);

%% Save MATLAB version
save('/path/to/whole_brain_TDM_mni_icbm152_nlin_asym_09b_new_script_meta_final.mat', 'TDM', '-v7.3');
fprintf('Done.\n');
