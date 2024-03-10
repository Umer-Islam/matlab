%function file is 'sort1.m'

clearvars;clc;

m = 2;
n = 6;
m1 = randn(m, n);

[unsortedVector, sortedVector] = func_A1Q4(m1);

disp('Original Matrix:');
disp(m1);

disp('Unsorted Column Vector:');
disp(unsortedVector);

disp('Sorted Column Vector:');
disp(sortedVector);
