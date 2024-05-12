%% int 
syms x 

f = sin(x);

y = int(f,x)

fplot(f)
%% 
clearvars,close all; clc, whos
syms x 
f = (x-1)^-1;

xi=0;
xf=2;
x = [xi,xf];
% int_f= int(f)
int_int = int (f, x, 0, 2, 'PrincipalValue', true)
fplot(f)


%%
clearvars,close all; clc, whos

