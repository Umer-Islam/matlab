%                                 Question#1
% 1. Write a MATLAB program which calculates the root of the equation c
% in the interval [0.5, 2], using False position method. Make sure to account for every
% possibility. Plot the given function and the calculated root.
clearvars;clc;who;

% Define the function
f = @(x) x^3 - 2*x^2 - 4*x + 8;

% Define the interval [a, b]
a = 0.5;
b = 2;

% Tolerance for stopping criteria
tol = 1e-6;

% Maximum number of iterations
maxIter = 1000;

% Initialize variables
iter = 0;
error = inf;
root = NaN;

% Check if the signs of f(a) and f(b) are opposite
if sign(f(a)) == sign(f(b))
    error('f(a) and f(b) must have opposite signs');
end

% Main loop of the False Position method
while error > tol && iter < maxIter
    % Increment iteration count
    iter = iter + 1;

    % Calculate the next approximation of the root using False Position method
    root_next = (a*f(b) - b*f(a)) / (f(b) - f(a));

    % Check if the root has been found
    if abs(f(root_next)) < tol
        root = root_next;
        break;
    end

    % Update interval [a, b]
    if sign(f(a)) * sign(f(root_next)) < 0
        b = root_next;
    else
        a = root_next;
    end

    % Calculate the error
    error = abs(root_next - root);

    % Update the root
    root = root_next;
end

% Display the root and number of iterations
disp(['Root of the equation: ', num2str(root)]);
disp(['Number of iterations: ', num2str(iter)]);

% Plot the function and the root
x = linspace(0.5, 2, 1000);
y = f(x);
plot(x, y, 'b-', 'LineWidth', 1.5);
hold on;
plot(root, f(root), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('f(x)');
title('False Position Method');
legend('Function', 'Root', 'Location', 'best');
grid on;
hold off;

%%                             Question#2
% 2. Write a MATLAB program which calculates the root of the function f (x) = x3−4x−9,
% using Secant method. Make sure to account for every possibility. Plot the given function
% and the calculated root.
clearvars;clc;who;
%%                             Question#3
% 
% 3. Write a MATLAB program which calculates the root of the equation x4 − x2 = 1,
% using the (fixed-point) Iteration method. Make sure to account for every possibility. Plot
% the given function and the calculated root.
clearvars;clc;who;
%%                             Question#4
% Write a MATLAB program which calculates the root of the function f (x, y) = x2 −y3
% using MATLAB’s fsolve function.
clearvars;clc;who;
%%                             Question#5
% Write a MATLAB program which calculates the roots of the polynomial x3 = 2 − x.
% Plot the given polynomial and the calculated roots.
clearvars;clc;who;