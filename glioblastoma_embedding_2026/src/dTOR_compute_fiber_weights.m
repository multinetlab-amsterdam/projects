% dTOR_compute_fiber_weights
% Script to seed the structural connectome with an ROI input, calculating a weight for each fiber streamline based on intersection with the ROI
% Adapted to run in sbatch from Suresh Joel, General Electric Global Research, April 2018

%%%

%%author__ = Mona Zimmermann
%%contact__ = m.l.m.zimmermann@amsterdamumc.nl
%%date__ = 2026/04/30
%%status__ = Finished 

%%%%%%%%%%%%%%%%%%%%
% Requirements     %
%%%%%%%%%%%%%%%%%%%%
%%% spm

function dTOR_compute_fiber_weights_main_for_sbatch(path, sub_id)


    % The following lines select the script input and require the user to specify this input as instructed
    load('/path/to/dTOR_fibers_vox_half_mm.mat'); % Provide path and filename for voxel-resampled connectome file (this determines the resolution in which the streamline sampling will be conducted)
    
    disp(path);
    disp(sub_id);
    
    vta_file = fullfile(path, sub_id, [sub_id,'_tumormask2MNI152NLin2009bAsym_res-05mm.nii.gz']); % Insert ROI file here (must be in unzipped nifti format and matching the resolution of the voxel-resampled file being used)
    disp(vta_file);
    %unzip .nii.gz file
    gunzip(vta_file);
    vta_file_unzipped = strrep(vta_file,'.gz','');
    
    v=spm_vol(vta_file_unzipped);
    im=spm_read_vols(v);
    im=im./max(im(:));
    im(im<0)=0;
    
    % The following lines initialize fiber weights to zero and then assign weights based on fiber-ROI intersection
    fibers_wt=zeros(length(fibers_vox),1);
    
    for i=1:length(fibers_vox)
        for j=1:size(fibers_vox{i},1)
            if (any(fibers_vox{i}(j,:)<=0) || fibers_vox{i}(j,1)>size(im,1) || fibers_vox{i}(j,2)>size(im,2) || fibers_vox{i}(j,3)>size(im,3))
                continue
            else
                fibers_wt(i)=fibers_wt(i) + im(fibers_vox{i}(j,1),fibers_vox{i}(j,2),fibers_vox{i}(j,3));
            end
        end
        if rem(i,10000)==0
            disp(['Fiber number:',num2str(i)])
        end
    end
    
    save([vta_file(1:end-7),'_fiber_wts.mat'],'fibers_wt');% was first end -4 but then the saved file does not have the correct name.

    fprintf('Yaaaaay subject is finished %s:\n', sub_id)


    
