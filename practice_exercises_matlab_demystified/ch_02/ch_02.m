%Ch2 Vectors and Matrices

% how to create a vector
a = [1,2,3,4,5] % row vector
% or 
b = [5;4;3;2;1] % column vector

% scalar multiple of a vector
3*a
5*b
% how to take transpose, just add variable_name'

b'

%vector addition and subtraction
a_scalar_multiple = 3*a;
sum_a = a_scalar_multiple + a;
subtract_a = a_scalar_multiple -a


% combining two vectors to make one larger one,concat

larger_a = [a,a_scalar_multiple]% notice the , or row vector
% and ; for column vector

% create vectors with uniform spacing between each element
% x= [initial_element: increment : final_element ]

x = [0:3:30]


% lets say you need 10 equal spaced values of exp(x)

x_data_points = 0:1:10 

exp_x_data_points = exp(x_data_points)




%%
clearvars,clc, close all, who
format short
% x = [100:-5:80];
x =[3;6;9];
magnitude_x  = sum(x.*x);
square_magnitude_x = sqrt(magnitude_x);
disp(square_magnitude_x)
%% magnitude of a complex vector
u =[1i;1+2i;4];
x= conj(u) .* u
y = sum(x)
sqrt (y)
%% buildin dot product command
clearvars,clc
u =[1i;1+2i;4];
a = [1;4;7;];
b = [2;-1;5];
c = dot(a,b)
u_dot = dot(u,u)
%% cross product
f_x = [1,2,3];
B = [2,3,4];
cross_product= cross(f_x,B)

%%
zeros(4,1)


%% how to find rank of a matrix
clearvars,clc
A = [0 1 0 2; 0 2 0 4];

r = rank(A);

fprintf('the rank is: %d \n',r)
%% how to take inverse of a matrix
clc,clearvars
A = [1,2,3; 4, 5,6;7,8,9];
inv(A);
%% inverse of a matrix whose determinant is zero
clc,clearvars,who
b = [1 2 3;0 2 2; 1 4 5];
b_inv = inv(b);
det_b = det(b);

fprintf('determinant of b %d', det_b)

fprintf('inverse of b is %d \n', b_inv)
%% 
clc, clearvars,who;

rand(4,3)