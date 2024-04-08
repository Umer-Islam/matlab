%% bisection method
clc,clearvars,close all
f = @ (x)  3*x.^2 -5;
a = 1;
b = 2;
min_step= 10^-4;
f_a = f(a);
f_b = f(b);

if f_a*f_b >0
    error("function has same sign at both endpoints")
end
[my_root, iterations] = bisection_function(f,a,b,min_step)
format_spec = 'the calcuated root is %8.6f\n\n';
fprintf(format_spec,my_root)

