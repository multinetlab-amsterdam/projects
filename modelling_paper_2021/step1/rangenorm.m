function A = rangenorm(matrix)
% matrix normaliseren door range
matrix2 = matrix;                                       
matrix2(1:size(matrix2,2)+1:(size(matrix2,2))^2)=[];    % diagonaal deleten
mx = max(max(matrix2));                                 % bereken maximale waarde
mn = min(min(matrix2));                                 % bereken minimale waarde
A = (matrix - mn)./(mx-mn);                             % normaliseren door range
A = A.*~eye(size(matrix,1),size(matrix,2));             % diagonaal op nul
end