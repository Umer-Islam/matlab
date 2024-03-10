% A1Q3, on a single graph make a plot of the function sin(x), cos(x) and
% tan(x) for 0 <= x <= 2pi. Give each curve color, legend and graph labels
clc;clearvars;
x_i = 0 ;
x_f = pi;
% disp(pi);
x = linspace(x_i,x_f);
sin_x = sin(x);
cos_x = cos(x);
tan_x = tan(x);
plot(x, sin_x,'r',x,cos_x,'bl',x, tan_x,'g');
legend('sin(x)', 'cos(x)', 'tan(x)');
xlabel('x');
ylabel('y');
title('graph of sin(x), cos(x), and tan(x) in range 0 to 2pi');
