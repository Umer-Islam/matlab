
clc; clearvars;
%create random matrix in range -1 to 1
matrixRandom = -1+2 * rand(5);
disp(matrixRandom);
% find where negative values are
matrixRandomNegative =matrixRandom(matrixRandom<0);
disp('negative values \n')
disp(matrixRandomNegative);

% find where positive values are
matrixRandomPositive = matrixRandom(matrixRandom>0);
fprintf('positive values \n');
disp(matrixRandomPositive);
