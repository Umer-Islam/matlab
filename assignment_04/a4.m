% assignment#4

%% Q1
clearvars,clc,close all;whos

% Define the IVP parameters
t0 = 0;      % Initial time
x0 = 1;      % Initial value of x at t = t0
t_final = 1; % The time at which we want to find the solution
h = 0.1;     % Step size

% Define the function f(t, x) = dx/dt
f = @(t, x) t ; % Example: replace this with any other differential equation

% Initialize variables
t = t0;
x = x0;

% Perform the RK4 iterations
while t < t_final
    if t + h > t_final
        h = t_final - t; % Adjust the last step to end exactly at t_final
    end

    % Calculate the k values
    k1 = h * f(t, x);
    k2 = h * f(t + h/2, x + k1/2);
    k3 = h * f(t + h/2, x + k2/2);
    k4 = h * f(t + h, x + k3);

    % Update x and t
    x = x + (k1 + 2*k2 + 2*k3 + k4) / 6;
    t = t + h;
end

% Display the RK4 result
disp(['RK4 x(', num2str(t_final), ') = ', num2str(x)])

%% Q2
clearvars,clc,close all; whos

% Define parameters
t0 = 0;         % Initial time
x0 = 0;         % Initial height
v0 = 0;         % Initial velocity
t_final = 10;   % Final time
h = 0.01;       % Time step size

% Initialize arrays to store the solutions
t_values = t0:h:t_final;
x_values = zeros(size(t_values));
v_values = zeros(size(t_values));

% Set initial conditions
x_values(1) = x0;
v_values(1) = v0;

% RK4 method to solve the system of ODEs
for i = 1:length(t_values) - 1
    t = t_values(i);
    x = x_values(i);
    v = v_values(i);

    m = 321 - 24 * t;
    if m <= 0
        break;
    end

    F = @(v) (5370 - 981 - (v^(3/2) / log(2 + v))) / m;

    k1_v = h * F(v);
    k1_x = h * v;

    k2_v = h * F(v + 0.5 * k1_v);
    k2_x = h * (v + 0.5 * k1_v);

    k3_v = h * F(v + 0.5 * k2_v);
    k3_x = h * (v + 0.5 * k2_v);

    k4_v = h * F(v + k3_v);
    k4_x = h * (v + k3_v);

    v_values(i + 1) = v + (k1_v + 2 * k2_v + 2 * k3_v + k4_v) / 6;
    x_values(i + 1) = x + (k1_x + 2 * k2_x + 2 * k3_x + k4_x) / 6;

    % Display the current step results on terminal
    fprintf('t = %.2f, x = %.4f, v = %.4f\n', t, x_values(i+1), v_values(i+1))

end

% Plot the results
figure;
subplot(2, 1, 1);
plot(t_values, x_values, 'b-');
xlabel('Time (s)');
ylabel('Height (x)');
title('Height vs. Time');
grid on;

subplot(2, 1, 2);
plot(t_values, v_values, 'r-');
xlabel('Time (s)');
ylabel('Velocity (v)');
title('Velocity vs. Time');
grid on;


%% Q3

clc,clearvars,close all; 
% Define the ODE as an anonymous function
ode = @(t, x) sin(x + t^2) ;

% Define the time span and initial condition
tspan = [-4, 4];
x0 = 0;

% Solve the ODE using ode45
[t, x] = ode45(ode, tspan, x0);

% Create a grid for the vector field
[x_grid, t_grid] = meshgrid(linspace(-4, 4, 20), linspace(-4, 4, 20));

% Compute the vector field
dxdt = sin(x_grid) + t_grid.^2;

% Normalize the vector field for better visualization
dxdt_norm = dxdt ./ sqrt(dxdt.^2 + 1);

% Display the values of t and x at each step
disp(' t       x');
disp('--------------');
for i = 1:length(t)
    fprintf(' %.2f    %.4f\n', t(i), x(i));
end

% Plot the solution and the vector field
figure;
hold on;

% Plot the vector field using quiver
quiver(t_grid, x_grid, ones(size(t_grid)), dxdt_norm, 'r');

% Plot the solution curve
plot(t, x, 'b', 'LineWidth', 2);

% Set plot labels and title
xlabel('t');
ylabel('x');
title('Solution and Vector Field of dx/dt = sin(x) + t^2');
xlim([-4 4]);
ylim([-4 4]);
grid on;
hold off;
%% Q4
clc,clearvars,close all; whos
% Define parameters
a = 1.5;    % Growth rate of prey in the absence of predators
c = 0.8;    % Death rate of predator in the absence of prey
b = 0.1;    % Effect of predators on prey growth rate
d = 0.1;    % Effect of prey on predator's growth rate
x0 = 30;    % Initial population of prey
y0 = 10;    % Initial population of predator

% Define the ODEs
ode = @(t, y) [a*y(1) - b*y(1)*y(2); -c*y(2) + d*y(1)*y(2)];

% Define the time span
tspan = [0, 20]; % Time span for simulation

% Solve the ODEs using ode45
[t, y] = ode45(ode, tspan, [x0; y0]);

% Plot populations as functions of time
figure;
subplot(2, 1, 1);
plot(t, y(:, 1), 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Prey Population');
title('Prey Population vs. Time');
grid on;

subplot(2, 1, 2);
plot(t, y(:, 2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('Predator Population');
title('Predator Population vs. Time');
grid on;

% Phase-plane plot
figure;
plot(y(:, 1), y(:, 2), 'm', 'LineWidth', 2);
xlabel('Prey Population');
ylabel('Predator Population');
title('Phase Plane Plot');
grid on;
