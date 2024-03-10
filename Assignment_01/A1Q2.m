
% generate random matrix in range;
matrixSize = [8,8];
randomMatrix = randi([-5,25],matrixSize);
disp(randomMatrix);

% check how many of them are zero;
fprintf('number of zero values: \n')
disp(numel(randomMatrix)- nnz(randomMatrix));
