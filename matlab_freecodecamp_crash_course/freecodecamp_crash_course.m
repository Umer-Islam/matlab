%% freecodecamp crash course

%% basic math opeartions
clc,clearvars;
x= 10 ;
y = 29;

sum = x +y

% type whos to see the names of variables in command windows
 % or who

 w = "whatever"; % single quotation creates string
 w_1 = 'whatever1'; % double quotation creates characters
 z = 4+x;
 %% Vectors and matrices
clc, clearvars

x = 1:10 % x = startingPoint :space/gap: end point , with unit space. if you check the 
% size using $whos it will be 1x10 a row with 10 columns

% transpose of a matrix , simply add variable_name' to take transpose

x'

% define arrays using linspace( buildin function) x = linspace(starting value, ending value, number of terms you want),
% automatically generates 100 values if not specified the third argument is
% not given i.e no. of terms
q = linspace(0,20,5)


% manually defining vectors

vec_1 = [12 ,12 ,12 ,12 ,12] % for row vector

vec_2 = [21;21;21;21;21] % for column vector


% A nxm matrix
A = [1,2;3,4;5,6] % a= [row;row;row]

% scalar operations on matrices
A + 2; % HERE SPACES ARE MUST ⚠

A * 3;

% A * A; % we should get a error becuase (3x2)(3x2) cannot be multiplied
A %(3x2)
A' %(2x3)

 A * A' 

square1 = linspace(0,100,101)
% square1 ^2 , where we will get error. for element operations we must use
% square1 .^2

square1.^2

% matrix of one
ones(3,1) %ones(rows,columns)
matrix_one = ones(4)
% similarly for zero's
matrix_zero= zeros(4)

% for identity matrix

matrix_identity = eye(4)

matrix_one + matrix_identity % no. of cols of first matrix must match the no. of rows of second

row_vector = [1,2.3,4.5,45] % only take out 4.5
row_vector(3)

index_matrix = [1,2,3,4; 4,5,6.221,333] % select 6
index_matrix(2,1)


% for the last value of matrix
index_matrix(end)
numel(index_matrix)

% to get the whole second row of index_matrix

index_matrix(2,1:2) % first we have selected the row then we select the the elements that we want

index_matrix(2,1:end)  % this is same as index_matrix(2,:)
index_matrix(2,:)
%% example problem
clc,clearvars,close all;
% A. waht is the max value of the following
% equation in the range 0<x<5
% y = -(x-3)^2 +10

x= linspace(0,5);

y = (-(x-3).^2) + 10;

plot(x,y,'*')
min(y)
max(y)

% use 'help min' for help
% use 'doc min' for documentation




