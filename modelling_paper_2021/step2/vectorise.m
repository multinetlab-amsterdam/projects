%%%%% vectorise adjacency matrix
function data = vectorise(matrix)
data =[];
% vector maken van x_matrix
for i=1:size(matrix,1)
    for j=1:i-1
        data=[data; matrix(i,j)];
    end
end
end