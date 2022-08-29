function cortical2surface(cortical, K, wsname, atlas)
% Project volumetric data parcellated according to the BNA atlas to the
% a fsaverage5 space. Required for Spin Test.
%
% FORMAT cortical2surface(cortical, wsname)
% cortical1    - the csv filename of cortical data to spin (format 210x2)
%                   the implementation is expecting to get a csv file.
% wsname       - the name of the directory to save the data
% K            - K is the number of vertices to project the value of the
%                 centroid. Recommened K=1 for use with SpinTest.m
%                 (centroid is a single vertex: then nearest neighnour to
%                 geometric centroid.) K = 0 inputs the value for the area
%                 to every vertex.
% atlas        - 'BNA' for brainnetome246 atlas, 'AAL' for automated
%                  anatomical labeling atlas. Default: 'BNA'.
% 
% OUT
% data.mat    - surface data. struct with 4 variables - surface[1,2]_[l,r]h
% for the 1st and 2nd datasets for the left and right hemispheres.
% exclude.mat - vertitices to exclude from downstream analyses. Vertices
% are excluded if these do not correspond to any specific BNA area (i.e.
% vertices are labeled as belonging to an unknown area) and if they belong
% to the medial wall.
% 
%
% EXAMPLE: cortical2surface('activity.csv','tumor_prob.csv', './tmp/', 'BNA')
% will project activity and tumor probability data on the surface and save
% the output files on ./tmp/
% 
%
%%author__ = Bernardo Maciel
%%contact__ = b.maciel@amsterdamumc.nl
%%date__ = 2021/11/22 % date script was created
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
% Expects ./atlas/ folder with freesurfer annotation and label files for
% both hemispheres for freesurfer and AAL atlases.
%%% Other m-files required:
% nearestneighbour.m (https://nl.mathworks.com/matlabcentral/fileexchange/12574-nearestneighbour-m)

% Input of atlases and definition of filenames
if atlas == 'AAL'
    ra = './atlas/rh.AAL_Atlas.annot';
    rl = './atlas/rh.AAL_cortex.label';
    la = './atlas/lh.AAL_Atlas.annot';
    ll = './atlas/lh.AAL_cortex.label';
else 
    ra = './atlas/rh.BN_Atlas.annot';
    rl = './atlas/rh.BN_cortex.label';
    la = './atlas/lh.BN_Atlas.annot';
    ll = './atlas/lh.BN_cortex.label';
end

% Read BNA annotation files
% *_labels holds the correspondence between the vertex and its
% annotation value. ex: vertex 1 left hemisphere = 8675760; vertex 1
% right hemisphere = 6388515.
[~, left_labels, ctl] = read_annotation(fullfile(la));
[~, right_labels,ctr] = read_annotation(fullfile(ra));

% Create wsname if it does not exist
wsname = string(wsname);
if ~exist(wsname, 'dir')
    mkdir(wsname)
end

% Load the cortical data files
data = load(cortical); data_spin1 = data(:,1); data_spin2 = data(:,2);

% Ensure uniqueness of every entry by manipulating the 6th decimal
for area = 1:length(data_spin1)
    data_spin1(area) = data_spin1(area) + area * 1e-6;
    data_spin2(area) = data_spin2(area) + area * 1e-6;
end

% Initialise variables to loop over
surface1_lh = nan(10242,1); surface2_lh = nan(10242,1);
surface1_rh = nan(10242,1); surface2_rh = nan(10242,1);

% ct? is a structure with several informations related to the key.
% ct?.table has the correspondence between BNA brain areas and
% annotation value. 8675760 = A6cdl_L (area 56); 6388515 = A6cdl_R
% (area 57)
bna_area_lh = nan(10242,1); bna_area_rh = nan(10242,1);
right_key = ctr.table(2:end,5);
left_key  = ctl.table(2:end,5);
x_left = false(10242,1); x_right = false(10242,1);

% for each vertex
for vertex=1:10242
    % get the brain area for each vertex of left and right hemispheres  
    % will get empty vector if not found or medial wall
    area_lh = find(left_key == left_labels(vertex));
    area_rh = find(right_key == right_labels(vertex));

    % populate surface vector with the value correspondent to the brain
    % area the vertex belongs to.
    if ~isempty(area_lh)
        bna_area_lh(vertex) = area_lh;
        surface1_lh(vertex) = data_spin1(area_lh);
        surface2_lh(vertex) = data_spin2(area_lh);
    else
        % if not found, mark vertex to exclude
        bna_area_lh(vertex) = NaN;
        x_left(vertex) = true;
    end

    % same for right
    if ~isempty(area_rh)
        bna_area_rh(vertex) = area_rh;
        surface1_rh(vertex) = data_spin1(area_rh);
        surface2_rh(vertex) = data_spin2(area_rh);
    else
        bna_area_rh(vertex) = NaN;
        x_left(vertex) = true;
    end    
end

if K > 0
    x_left = true(10242,1); x_right = true(10242,1);
    % Import the geometric data from the location of the vertices of the
    % left hemisphere
    fileID = fopen(ll, 'r');
    coords_left  = cell2mat(textscan(fileID, '%f %f %f %f %f', 'HeaderLines', 2));
    fclose(fileID);
    % Same for the right hemisphere
    fileID = fopen(rl,'r');
    coords_right = cell2mat(textscan(fileID, '%f %f %f %f %f', 'HeaderLines', 2));
    fclose(fileID);

    for area=1:length(data_spin1)
        % Find all of the certices corresponding to a BNA area
        left_hits = find(bna_area_lh == area);
        % If the area exists is big enough on the right side
        if ~isempty(left_hits) && length(left_hits) > K
            % Find all of the vertices coordinates
            subset_lh = coords_left(ismember(coords_left(:, 1), left_hits), :);
            % Calculate the (x,y,z) coordinates of the real centroid of a
            % region
            xx = mean(subset_lh(:, 2));
            yy = mean(subset_lh(:, 3));
            zz = mean(subset_lh(:, 4));
            % Find the K-nearest-neighbours to the centroid to atribute the
            % value of the data.
            nn = nearestneighbour([xx yy zz]', subset_lh(:, 2:4)', 'NumberOfNeighbours', K);
            centroid_lh = subset_lh(nn, 1);
            % Not exclude your centroids
            x_left(centroid_lh) = false;
        end
        % Same for right side
        right_hits = find(bna_area_rh == area);
        if ~isempty(right_hits) && length(right_hits) > K
            subset_rh = coords_right(ismember(coords_right(:, 1), right_hits), :);

            xx = mean(subset_rh(:, 2));
            yy = mean(subset_rh(:, 3));
            zz = mean(subset_rh(:, 4));

            nn = nearestneighbour([xx yy zz]', subset_rh(:, 2:4)', 'NumberOfNeighbours', K);
            centroid_rh = subset_rh(nn, 1);
            x_right(centroid_rh) = false;
        end
    end
end

save(strcat(wsname, 'data.mat')   , 'surface1_rh', 'surface1_lh', 'surface2_rh', 'surface2_lh')
save(strcat(wsname, 'exclude.mat'), 'x_left'     , 'x_right')
