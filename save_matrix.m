%% Save_matrix.m
% save a matrix in csv format
% Dos Reis F.
% 02.2021
function save_matrix(file,matrix)
    writematrix(size(matrix),file);
    writematrix(matrix,file,'Delimiter',',','WriteMode','append');
end