clc,clearvars,whos, close all 
% define the function
f = @(x) x^3 -2*sin(x);
% define the interval
a = 0.5;
b= 2; 
% define the itterations
n = 29;
e = 0.00001;
% write conditions
if f(a)*f(b)<0
	for i = i :n
		c = (a+b)/2;
        fprintf('root is : %.5f ,  no. of itteration:  %d\n',c,i)
        if abs(a-c)<e || abs(b-c)<e
            break
        end
		if f(a)*f(c)<0
		b =c;
		else 
		a =c;
		end
	end
else
disp('No root lies in the interval')
end
%plot function and its root in the given interval
hold on 
xlabel("x",Interpreter="latex",FontSize=15)
ylabel("f(x)",Interpreter="latex",FontSize=15)
legend('Interpreter','latex',FontSize=15) 
box off
fplot(f,lineWidth=1.5)
plot(c,0,'*r',LineWidth=2)
grid on 
hold off