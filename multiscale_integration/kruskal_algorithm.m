function [T,w,PP] = kruskal_algorithm(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kruskal's algorithm to construct an MST
% Main part of the code
% X is weighted adjacency matrix, assuming that X is a connected matrix!
% output T --> MST adjacency matrix

trilmatrix = tril(X,-1);
[r, c, v] = find(trilmatrix);
PV = [r c v];
n = size(X,1);
row = size(PV,1);

% sort PV by ascending weights order. 
[PPV, in] = sort(PV(:,3));
rowe = PV(in,1);
col = PV(in,2);
PB = [rowe col PPV];
korif = zeros(1,n);
PP = flipud(PB);

T = zeros(n);
for i = 1 : row
% control if we insert edge[i,j] in the graphic. Then the graphic has
% circle
    akmi = PP(i,[1 2]);
    [korif,c] = iscycle(korif,akmi);
    if c == 1
       PP(i,:) = [0 0 0];
   end
end
% Calculate Minimum spanning tree's weight
w = sum(PP(:,3)');
% Create minimum spanning tree's adjacency matrix
for i = 1 : row
    if PP(i,[1 2]) ~= [0 0]
        T(PP(i,1),PP(i,2)) = PP(i,3);
        T(PP(i,2),PP(i,1)) = PP(i,3);
    end
end
    T = T>0;
    T = double(T);
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
g=max(korif)+1;
c=0;
n=length(korif);
if korif(akmi(1))==0 & korif(akmi(2))==0
    korif(akmi(1))=g;
    korif(akmi(2))=g;
elseif korif(akmi(1))==0
    korif(akmi(1))=korif(akmi(2));
elseif korif(akmi(2))==0
    korif(akmi(2))=korif(akmi(1));
elseif korif(akmi(1))==korif(akmi(2))
    c=1;
    return
else
    m=max(korif(akmi(1)),korif(akmi(2)));
    for i=1:n
        if korif(i)==m
           korif(i)=min(korif(akmi(1)),korif(akmi(2)));
       end
   end
end
end