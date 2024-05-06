% A3
%% Question#01
%% Question#02
%% Question#03
%% Question#04
%% Question#05
%% Question#06  Evaluate the following double integral using MATLAB’s int function.
clearvars,clc,close all;whos

f = @(x,y) (x.^2-10*x.*y-y.^2)^4;
x= linspace(0,1);
y= linspace(0,1);
%intergra2  ...Q = integral2(FUN,XMIN,XMAX,YMIN,YMAX)...
% x_lim = [0,1];
% y_lim = [0,1];
% integral_2_buildin = integral2(f,0,1,0,1);
% plot(f,x,y)


%dblquad ...Q = dblquad(integrnd, pi, 2*pi, 0, pi)...
double_quad = dblquad(f,1,2,1,2)

% trapz


%% test
clc,close all; clearvars, whos
% Define the function f
f = @(x,y) (x.^2-10.*x.*y-y.^2).^4;

x_lim= [0,1];
y_lim = [0,1];

%intergral2
tic;
integral_2 = integral2(f,x_lim(1),x_lim(2),y_lim(1),y_lim(2));
time_integral2 = toc

%dlbquad
tic;
dbl_quad= dblquad(f,x_lim(1),x_lim(2),y_lim(1),y_lim(2));
time_dblquad= toc;
%trapz
