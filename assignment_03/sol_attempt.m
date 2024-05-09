% A3
%% Question#01
clearvars,clc,close all;whos
x = [1, 3/2, 0, 2];
y = [3, 13/4, 3, 5/3];
degree = length(x) - 1;
p = polyfit(x, y, degree);

fprintf("Coefficients of the interpolating polynomial: ")
disp(p)
%to turn cofficient into equation
syms x

polynomial_expr = poly2sym(p, x);
fprintf("Equation for the data points given: ")
disp(polynomial_expr)
first_derivative = diff(polynomial_expr, x);
fprintf("First derivative: ")
disp(first_derivative)
second_derivative = diff(first_derivative, x);
fprintf("second derivative: ")
disp(second_derivative)


hold on
fplot(polynomial_expr,'-r')
fplot(first_derivative,'-g')
fplot(second_derivative,'-b')
legend('Polynomial', 'First Derivative', 'Second Derivative','interpreter','latex');

grid on
hold off

%% Question#02
% trapizod rule
clearvars, clc, close all; whos

% f = @(x) (x^-3);
f = @(x) x+0.5;
a = 0;
b=2;
n = 2;
h = (b-a)/n;
s = 0.5 * (f(a)+ f(b));

for i = 1:n-1
    s = s + f(a + i*h);
end

I = h * s;
fprintf("Approximate value of integral is:")
disp(I)
%% Question#03
% simpson rule
clc,close all; clearvars,whos
f = @(x) x+0.5;
% f = @(x) x^-3;
a = 0 ;
b = 2;
n = 2;
if mod(n,2)~=0
    error("Cannot use this method number of intervals are odd")
end
h = (b-a)/b;
s = f(a)+f(b);
for i = 1:2:n-1
    s = s+4*f(a+i*h);
end
for i = 2:2:n-2
    s =s + 2*f(a+i*h);
end
I = h/3 *s;
fprintf("Approximate value of integral is:")
disp(I)

%% Question#04
clc,clearvars, close all; whos
% range
x_range = 1:0.1:3;  
y_range = 0:0.1:3;  
z_range = 0:0.1:1;

% meshgrid
[X, Y, Z] = ndgrid(x_range, y_range, z_range);

% function/integrand
f = @(x, y, z) x.^2 + y.^2 + z.^2;

F = f(X, Y, Z);

I = trapz(z_range, trapz(y_range, trapz(x_range, F, 1), 2), 3);

fprintf("Calculated triple integral value:  %.5f\n", I);

%% Question#05
clc,clearvars,close all; whos
%trapz
tic;
x = 0:.1:1; 
y = 0:.1:1; 
[X,Y] = meshgrid(x,y);
F =  X.^2- 10 *X.*Y + Y.^2;
I = trapz(y,trapz(x,F,2));
fprintf('trapz value: %.6f\n', I);

trapz_time= toc;

%integral2
f = @(x,y) x.^2- 10 *x.*y + y.^2;
tic;
integral_2 = integral2(f,0,1,0,1);
fprintf('integral_2 value: %.6f\n', integral_2);

integral_2_time = toc;
%dblquad
tic;
double_quad = dblquad(f,0,1,0,1);
fprintf('double_quad value: %.6f\n', double_quad);

double_quad_time= toc;
%fastest
lowest_value = min([trapz_time, double_quad_time, integral_2_time]);
if lowest_value== lowest_value(1)
    fprintf('trapz produced the lowest value.\n');
elseif lowest_value== lowest_value(2)
    fprintf('dblquad produced the lowest value.\n');
elseif lowest_value == lowest_value(3)
    fprintf('intagral2 produced the lowest value.\n');
end


%% Question#06  
% Evaluate the following double integral using MATLAB’s int function.
clc,clearvars,close all; whos

syms x t

fx = x^2;
ft =  3 * log(t);

x_start = 1;
x_end = 2;
t_start = 0;
t_end = 1;

fx_int = int(fx, x, x_start, x_end);

ft_int = int(ft, t, t_start, t_end);
II = fx_int+ft_int;

fprintf('The double integral of f = x^2 + 3 * log(t) over the range x = [0,1] and t = [1,2] is: %s\n', II);
