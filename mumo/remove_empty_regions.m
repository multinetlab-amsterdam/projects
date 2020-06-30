function [clean_data_1, deleted_regions_1, clean_data_2, deleted_regions_2] = remove_empty_regions(raw_data_1, raw_data_2)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % 02.JUNE.2020
   % L C BREEDT
   % 
   %
   % This script checks for and removes empty regions in a dataset when
   % only one input argument is provided; or two datasets (e.g. fMRI & DWI)
   % when two input arguments are provided.
   %
   %
   % Input:
   %     A dataset 1 and an optional dataset 2 that need(s) to be checked
   %     for empty regions.
   %
   % Output:
   %    A clean dataset 1 and a clean dataset 2, as well as a list for each
   %    dataset that specifies which regions have been removed.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize variables
    nrois = length(raw_data_1(1,1,:));
    nsubs_1 = length(raw_data_1(:,1,1));
    deleted_regions_1 = [];
    
    % create atlas reference - atlas INDEX is "new", updated atlas
    % region; atlas VALUE is original, actual atlas region
    atlas_ref = 1:nrois;
    
    % implement exception in case only one dataset is given as input
    % argument
    if nargin < 2
        nsubs_2 = length(raw_data_1(:,1,1));
        raw_data_2 = ones(nsubs_2, nrois, nrois);
    end
    
    % initialize additional variables in case two datasets are given as
    % input
     nsubs_2 = length(raw_data_2(:,1,1));
     deleted_regions_2 = [];
     
    % can we use a "while" loop? so first initialize two variables (one for
    % each dataset) that are non-empty as long as there are still 0s in one
    % of the datasets; then, WHILE variable 1 OR variable 2 are non-empty,
    % the function will run. it will start with dataset 1, removing regions
    % IF variable 1 is non-empty - from both datasets - and updating
    % variable 1 and 2. IF variable 1 is empty, it will remove regions IF
    % variable 2 is non-empty - from both datasets - and update variable 1
    % and 2. when both variables are empty, the removal is done and the
    % script will output two clean datasets and a list of deleted regions
    
    %% script body
    

    
    % initialize conditionals
    conditional_1 = [];
    for s = 1:nsubs_1
        if ~isempty(find(all(raw_data_1(s,:,:)==0)))
            currsub_empty_rois = find(all(raw_data_1(s,:,:)==0));
            if isempty(conditional_1)
                conditional_1 = currsub_empty_rois; % initialize vector containing all empty rois
            end
            conditional_1 = vertcat(conditional_1, setdiff(currsub_empty_rois, conditional_1));
        end
    end
    
    conditional_2 = [];
    for s = 1:nsubs_2
        if ~isempty(find(all(raw_data_2(s,:,:)==0)))
            currsub_empty_rois = find(all(raw_data_2(s,:,:)==0));
            if isempty(conditional_2)
                conditional_2 = currsub_empty_rois; % initialize vector containing all empty rois
            end
            conditional_2 = vertcat(conditional_2, setdiff(currsub_empty_rois, conditional_2));
        end
    end
    
    
    % while-loop that runs until there are no more zeroes in either dataset
    while ~isempty(conditional_1) | ~isempty(conditional_2)
        if ~isempty(conditional_1)
            empty_rois = [];
            for s = 1:nsubs_1
                if ~isempty(find(all(raw_data_1(s,:,:)==0)))
                    currsub_empty_rois = find(all(raw_data_1(s,:,:)==0));
                    if isempty(empty_rois)
                        empty_rois = currsub_empty_rois; % initialize vector containing all empty rois
                    end
                    empty_rois = vertcat(empty_rois, setdiff(currsub_empty_rois, empty_rois));
                end
            end
            empty_rois = sort(empty_rois);
    
            % delete dataset 1 empty regions in both datasets
            raw_data_1(:,:,empty_rois) = [];
            raw_data_1(:,empty_rois,:) = [];
            
            raw_data_2(:,:,empty_rois) = [];
            raw_data_2(:,empty_rois,:) = [];
    
            % save empty rois in deleted_regions for output - based on actual atlas
            % regions, so based on atlas_ref!
            actual_regions = atlas_ref(empty_rois)';
            deleted_regions_1 = vertcat(deleted_regions_1, actual_regions);
    
            % update atlas reference to reflect region removal
            atlas_ref(empty_rois) = [];
            
        else % if conditional_1 is no longer empty, conditional 2 still is, so "else" (since whole loop hangs on conditional_1 OR 2 being empty)
            empty_rois = [];
            for s = 1:nsubs_2
                if ~isempty(find(all(raw_data_2(s,:,:)==0)))
                    currsub_empty_rois = find(all(raw_data_2(s,:,:)==0));
                    if isempty(empty_rois)
                        empty_rois = currsub_empty_rois; % initialize vector containing all empty rois
                    end
                    empty_rois = vertcat(empty_rois, setdiff(currsub_empty_rois, empty_rois));
                end
            end
            empty_rois = sort(empty_rois);
    
            % delete dataset 1 empty regions in both datasets
            raw_data_1(:,:,empty_rois) = [];
            raw_data_1(:,empty_rois,:) = [];
            
            raw_data_2(:,:,empty_rois) = [];
            raw_data_2(:,empty_rois,:) = [];
    
            % save empty rois in deleted_regions for output - based on actual atlas
            % regions, so based on atlas_ref!
            actual_regions = atlas_ref(empty_rois)';
            deleted_regions_2 = vertcat(deleted_regions_2, actual_regions);
    
            % update atlas reference to reflect region removal
            atlas_ref(empty_rois) = [];
        end
        
        % update conditionals
        conditional_1 = [];
        for s = 1:nsubs_1
            if ~isempty(find(all(raw_data_1(s,:,:)==0)))
                currsub_empty_rois = find(all(raw_data_1(s,:,:)==0));
                if isempty(conditional_1)
                    conditional_1 = currsub_empty_rois; % initialize vector containing all empty rois
                end
                conditional_1 = vertcat(conditional_1, setdiff(currsub_empty_rois, conditional_1));
            end
        end
    
        conditional_2 = [];
        for s = 1:nsubs_2
            if ~isempty(find(all(raw_data_2(s,:,:)==0)))
                currsub_empty_rois = find(all(raw_data_2(s,:,:)==0));
                if isempty(conditional_2)
                    conditional_2 = currsub_empty_rois; % initialize vector containing all empty rois
                end
                conditional_2 = vertcat(conditional_2, setdiff(currsub_empty_rois, conditional_2));
            end
        end
        
    end
    
    % create export variables
    clean_data_1 = raw_data_1;
    clean_data_2 = raw_data_2;
    
    deleted_regions_1 = sort(deleted_regions_1);
    deleted_regions_2 = sort(deleted_regions_2);
    
end