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




