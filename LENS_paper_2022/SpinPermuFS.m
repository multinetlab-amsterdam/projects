function SpinPermuFS(datal,datar,permno,wsname)
% Compute designated # of permutations/spins of the input surface data
% in FreeSurfer fsaverage5.
% FORMAT SpinPermuFS(readleft,readright,permno)
% readleft     - left surface data to spin 
% readright    - right surface data to spin
% permno       - the number of permutations
% wsname       - the name of a workspace file including all spun data to be saved
% Example   SpinPermuFS('../data/depressionFSdataL.csv','../data/depressionFSdataR.csv',100,'../data/rotationFS.mat')
% will spin prebuilt data, neurosynth map associated with 'depression', 100
% times, and save the workspace file of all spun data in ../data/rotationFS.mat
% Aaron Alexander-Bloch & Siyuan Liu 
% SpinPermuFS.m, 2018-04-22
% The implementation of generating random rotations originally described in our paper — 
% rotating the coordinates of vertices at angles uniformly chosen between zero and 360 degrees
% about each of the x (left-right), y (anterior-posterior) and z (superior-inferior) axes —
% introduces a preference towards oversampling certain rotations. 
% Thus, we modified the code to incorporate an approach, Lefèvre et al. (2018), 
% that samples uniformly from the space of possible rotations. The updated
% uniform sampling prodcedure does not require AxelRot.m anymore.
% Updated on 2018-07-18
% Update 07/31/2020 (SMW): will automatically remove medial wall for
% fsaverage5. may need to change if not fsaverage5 (10242 vertices per
% hemisphere)
% Update 19-11-2021 Bernardo Maciel: the script was adapted to have as
% inputs the surface data with NaN on the vertices to remove.


%Set up paths
fshome = getenv('FREESURFER_HOME');
fsmatlab = sprintf('%s/matlab',fshome);
path(path,fsmatlab);


% extract the corresponding sphere surface coordinates for rotation
[verticesl, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage5/surf/lh.sphere'));
[verticesr, ~] = freesurfer_read_surf(fullfile(fshome,'subjects/fsaverage5/surf/rh.sphere'));

% Use rng to initialize the random generator for reproducible results.
rng(0);

%initialize variables to save rotation
bigrotl= nan(permno, 10242);
bigrotr= nan(permno, 10242);

I1 = eye(3,3);
I1(1,1)=-1;
bl=verticesl;
br=verticesr;

%permutation starts
for j=1:permno
    if rem(j, 1000)==0
        j
    end
    %the updated uniform sampling procedure
    A = normrnd(0,1,3,3);
    [TL, temp] = qr(A);
    TL = TL * diag(sign(diag(temp)));
    if(det(TL)<0)
        TL(:,1) = -TL(:,1);
    end
    %reflect across the Y-Z plane for right hemisphere
    TR = I1 * TL * I1;
    bl = bl*TL;
    br = br*TR;    
    
    %Find the pair of matched vertices with the min distance and reassign
    %values to the rotated surface.
    Il = nearestneighbour(verticesl', bl'); % added 2019-06-18 see home page
    Ir = nearestneighbour(verticesr', br'); % added 2019-06-18 see home page

    %save rotated data
    bigrotl(j, :) = datal(Il);
    bigrotr(j, :) = datar(Ir);
end
save(wsname,'bigrotl','bigrotr')
%save bigrotl and bigrotr in a workspace file for the null distribution
