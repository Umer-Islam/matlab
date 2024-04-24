%% plot a simple trig function
clearvars, clc, close all;

% define the interval on which you want to see the function plotted,make
% sure that data points are defined before the function

x = [-6:0.1:6];
f = sin(x);

plot(x,f)
xlabel("x")
ylabel("y")

%% another plot
clearvars,clc,close all;
t = 0:0.02:4 ;
f = exp(-2*t).*sin(t);

g = exp(-t).*sin(t);
plot(t,f,'-',t,g,'--')
% fplot('exp(-2.*t).*sin(t)',[0,4])
title("damped spring force")
legend('f','g')
grid on
% axis equal
axis auto

%% setting axis in graph
clearvars,clc,close all;
x = [0:0.1:5];
y = sin(2*x + 3);
plot(x,y)
axis([-10,10,-20,20])

%% plotting
clearvars, clc,close all;
y = exp((-1.5)*x).*sin(5*x+3);
plot(x,y)
axis([0, 5, -1, 1])
%% subPlot
clearvars , clc, close all
x = [0:0.01:5];
x1 = 5:0.01:10;
y = exp(-1.2*x).*sin(20*x);
subplot(1,2,1)
plot(x,y)
%% polar plot
clearvars, clc, close all
theta = 0: pi/90: 6*pi;
a=2;
r = 1 + a * cos(theta);
polarplot(theta,r,'r ')
%%  loglog plot
clearvars, clc, close all
 x = [0:0.1:20];
y = exp(-10*x.^2);
% plot(x,y)
loglog(x,y)
legend('plot','loglog')
%% graph of 2 types of data points
clearvars, clc, close all
student = [1:6];
pointsInSubject = [45,99,65,34,23,65];
plot(student,pointsInSubject)
bar(student,pointsInSubject,'cyan')
%% stem plot
clearvars ,clc ,close all;
t = [0:5:200];
f = exp(-0.01*t).*sin(t/4);
stem(t,f,'dr')