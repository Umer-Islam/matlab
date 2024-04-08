% %  matrix operations
clc,clearvars,close all
a = magic(4)  %4x4
b = magic(6);
random_matrix = rand(4,6); %4x6
% c = [a,b] % error due to differnt dimentions

% d = [a;b] %error due to inconsistent dimentions

e = [a,random_matrix];  %row wise concat possible 
% due to same no. of rows

% f = [a;random_matrix] % col wise concat NOT
% possible due to different no. of cols

g = a(:,3)
h = a(1:3,3)

% manual matrix defining
a1 = [1:3] % row wise

b1 = 5:7;

a2 = dot(a1,b1);
a3 = dot(a1,b1');

a4 = cross(a1,b1)

a5 = kron(a1,b1)
a6 = length(a1);
a7 = max(b1)
a8 = min(b1)

A = [1,2,3,4;5,6,7,8;9,10,11,12]

size(A)
length(A)
% A(1,1) =20
A
A(:,2)

sum(A) % adds all elements of a column
A
max(A)

%%
clc,clearvars,close all
% create 5x5 magic matrix, verify that sum of
% intagers in each row, column and diagonal is equal
A = magic(5)
sum(A)
sum(A,2)
% both do the same thing 
% sum(diag(A)) == trace(A)
sum(diag(A))
trace(A)
flip_A=flip(A)
sum(diag(flip_A))

A(3,:)
A
max(max(A))

A
sort(A);% sorts each column lowest to highest
sort(A,2)% sorts each row lowest to highest

sort(A,2,"descend")
det(A)

%% for Loop
% for loop
clc,clearvars,close all
my_sum = 0;
for i = 1:3
    A(i+1) = 1
end
sum(A)

%% if else loop
a =6; 
if  a >1 && a<10 
    disp('true')
% elseif a == 7
%     disp('equal')
% else 
%     disp('whatever')
end


%% while loop
clc,clearvars,close all
i = 0;
while i<10
disp('1')
i = i+1;
pause
end
disp('0')

%% 2d plotting 
clc,clearvars,close all
x_i = 0;
x_f = 4* pi;
x = linspace(x_i,x_f)
sin = sin(x);
cos = cos(x);
plot(x,sin,x,cos,'+')
xlim([x_i,x_f])
ylim ([-2,2])
xlabel('$$x$$','Interpreter','latex')
ylabel('$$y$$','Interpreter','latex')
legend({'$$\sin(x)$$','$$\cos(x)$$'},'Interpreter','latex','fontsize',9)
legend boxoff
%% 3d plot
clc,clearvars,close all
x = linspace(0,10*pi,200)
sin = sin(x);
cos = cos(x);
plot3(sin,cos,x,'lineWidth',2)
xlabel('sinx')
ylabel('cosx')
zlabel('x')
grid on
title("helix")

%% surface plot exercise
 clc,clearvars,close all

x = linspace(-5,5);
y = linspace(-5,5);
[X,Y] = meshgrid(x,y)
z = Y.*sin(X)- X*cos(Y)

surf(X,Y,z)



