function [T,w,PP] = kruskal_algorithm(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kruskal's algorithm to construct an MST
% Main part of the code
% X is weighted adjacency matrix, assuming that X is a connected matrix!
% output T --> MST adjacency matrix

trilmatrix = tril(X,-1); % lower triangle of matrix X
[r, c, v] = find(trilmatrix); % find rowindex, columnindex, and value at that index of all values in matrix
PV = [r c v]; % place rowindex, columnindex, and value at that index in matrix PV
n = size(X,1); % size of matrix, so how many nodes
row = size(PV,1); % size of PV, so how many edges (or how many rows)

% sort PV by ascending weights order. 
[PPV, in] = sort(PV(:,3)); % sort weights (edges) in PV in ascending order. PPV = PV sorted; in = original indices (so index of relevant weight in PV)
rowe = PV(in,1); % rowindex of sorted weights (so rowe(1) = rowindex of weakest weight as sorted in PPV)
col = PV(in,2); % columnindex of sorted weights
PB = [rowe col PPV]; % rowindex (so node), columnindex (so node), weight in ascending order (weakest first)
korif = zeros(1,n); % pre-allocate empty list for later, see iscycle
PP = flipud(PB); % flip PB upside down - so strongest weight is now at top of list and vice versa

% determine indices of duplicate values
A = PP(:,3);
[~, unique_idx] = unique(A); % returns list of all unique values in A
duplicate_id = setdiff(1:numel(A), unique_idx); % returns indices of duplicate values in A

T = zeros(n); % empty matrix (pre-allocated MST matrix)
for i = 1 : row % for all edges
% control if we insert edge[i,j] in the graphic. Then the graphic has
% circle
    akmi = PP(i,[1 2]); % if i=1, akmi is row- and columnindex of first entry in PP, so of strongest weight
    [korif,c] = iscycle(korif,akmi); % see iscycle below. korif is vector indicating for every node if it has links, c is variable indicating circle yes or no
    if c == 1 % if circle will be created by placing link, then: 
       PP(i,:) = [0 0 0]; % set rowindex, columnindex, and weight to 0 -> so no link will be formed, 0 in MST!
   end
end
% Calculate Minimum spanning tree's weight
w = sum(PP(:,3)'); % not further used, but this is the sum of all weights remaining after setting links forming circles to 0
% Create minimum spanning tree's adjacency matrix
for i = 1 : row % for all edges
    if PP(i,[1 2]) ~= [0 0] % if edge is not 0 (row- & columnindex not 0), then:
        T(PP(i,1),PP(i,2)) = PP(i,3); % set edge in T to edgeweight (Bxy = weight)
        T(PP(i,2),PP(i,1)) = PP(i,3); % set edge in T to edgeweight (Byx = weight)
    end
end
    T = T>0; % set edges to 1 (all values larger than 0)
    T = double(T); % double precision
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iscycle
%
% input: korif = set of vertices in the graph
%       akmi = edge we insert in graph
% output: korif = The "new: set of vertices
%        c = 1 if we have circle, else c = 0
%
% N.Cheilakos,2006
function [korif,c]=iscycle(korif,akmi)
g=max(korif)+1; % highest value of korif (so amount of links) + 1 (1 link added)
c=0; % c is boolean used to indicate cycle yes or no
n=length(korif); % amount of nodes
if korif(akmi(1))==0 & korif(akmi(2))==0 % if node akmi(1) and node akmi(2) both have no links yet, then:
    korif(akmi(1))=g; % node akmi(1) is set to g, so a link is placed
    korif(akmi(2))=g; % node akmi(2) is set to g, so a link is placed
elseif korif(akmi(1))==0 % if only node akmi(1) has no links (but akmi(2) does), then:
    korif(akmi(1))=korif(akmi(2)); % node akmi(1) is set to same value as akmi(2), so link is placed
elseif korif(akmi(2))==0 % if only node akmi(2) nas no links (but akmi(1) does), then:
    korif(akmi(2))=korif(akmi(1)); % node akmi(2) is set to same value as akmi(1), so link is placed
elseif korif(akmi(1))==korif(akmi(2)) % if akmi(1) and akmi(2) have same amount of links, then:
    c=1; % set c to 1, so placing link would form cycle?
    return
else % if none of the above applies (so both nodes have links but not equal)
    m=max(korif(akmi(1)),korif(akmi(2))); % m= most links either of the nodes has
    for i=1:n % for all nodes
        if korif(i)==m % if node i has same amount of links as max of akmi(1) and akmi(2), then:
           korif(i)=min(korif(akmi(1)),korif(akmi(2)));  % node i is set to lowest amount of links
       end
   end
end
end