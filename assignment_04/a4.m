% assignment#4

%% Q1
clearvars,clc,close all;whos

% Define the IVP parameters
t0 = 0;     
x0 = 1;      
t_final = 1; 
h = 0.5;    


f = @(t, x) t ; 


t = t0;
x = x0;


while t < t_final
    if t + h > t_final
        h = t_final - t; 
    end

   
    k1 = h * f(t, x);
    k2 = h * f(t + h/2, x + k1/2);
    k3 = h * f(t + h/2, x + k2/2);
    k4 = h * f(t + h, x + k3);


    x = x + (k1 + 2*k2 + 2*k3 + k4) / 6;
    t = t + h;
end

% Display the RK4 result
disp([' value of  x(', num2str(t_final), ') = ', num2str(x)])

%% Q2
clearvars,clc,close all; whos


t0 = 0;         
x0 = 0;         
v0 = 0;         
t_final = 10;   
h = 0.5;      


t_values = t0:h:t_final;
x_values = zeros(size(t_values));
v_values = zeros(size(t_values));

% ic
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

   
    fprintf('t = %.2f, x = %.4f, v = %.4f\n', t, x_values(i+1), v_values(i+1))

end


figure
subplot(2, 1, 1)
plot(t_values, x_values, 'b-')
xlabel('Time (s)','interpreter','latex')
ylabel('Height (x)','interpreter','latex')
title('Height vs. Time','interpreter','latex')
grid on

subplot(2, 1, 2)
plot(t_values, v_values, 'r-')
xlabel('Time (s)','interpreter','latex')
ylabel('Velocity (v)','interpreter','latex')
title('Velocity vs. Time','interpreter','latex')
grid on


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
dxdt = sin(x_grid + t_grid.^2);

% Normalize the vector field for better visualization
dxdt_norm = dxdt ./ sqrt(dxdt.^2 + 1);

disp(' t    |   x');
disp('--------------');
for i = 1:length(t)
    fprintf(' %.2f  |  %.4f\n', t(i), x(i))
end

figure
hold on

quiver(t_grid, x_grid, ones(size(t_grid)), dxdt_norm, 'r')

% Plot the solution curve
plot(t, x, 'b', 'LineWidth', 2)

% Set plot labels and title
xlabel('t','interpreter','latex')
ylabel('x','interpreter','latex')
title('Solution and Vector Field of dx/dt = sin(x + t^2)','interpreter','latex')
xlim([-4 4])
ylim([-4 4])
grid on
hold off
%% Q4
clc,clearvars,close all; whos
a = 1.5;    % absence of predators
c = 0.8;    % Death rate of predator in the absence of prey
b = 0.1;    % Effect of predators on prey growth rate
d = 0.1;    % Effect of prey on predator's growth rate
x0 = 30;    % Initial prey pop
y0 = 10;    % Initial predator pop

ode = @(t, y) [a*y(1) - b*y(1)*y(2); -c*y(2) + d*y(1)*y(2)];

tspan = [0, 20]; 

[t, y] = ode45(ode, tspan, [x0; y0]);

% Plot
figure
subplot(2, 1, 1)
plot(t, y(:, 1), 'b-')
xlabel('Time','interpreter','latex')
ylabel('Prey Population','interpreter','latex')
title('Prey Population vs. Time','interpreter','latex')
grid on

subplot(2, 1, 2)
plot(t, y(:, 2), 'r-')
xlabel('Time','interpreter','latex')
ylabel('Predator Population','interpreter','latex')
title('Predator Population vs. Time','interpreter','latex')
grid on

% Phase-plane plot
figure
plot(y(:, 1), y(:, 2), 'm', 'LineWidth', 2)
xlabel('Prey Population','interpreter','latex')
ylabel('Predator Population','interpreter','latex')
title('Phase Plane Plot','interpreter','latex')
grid on
%% Q5
clc,close all, clearvars, whos
syms x(t) omega


ode = diff(x, t, 2) + omega^2 * x == 0;
sol = dsolve(ode);
disp('General Solution:')
disp(sol)
x0 = sym('x0');
v0 = sym('v0');
Dx = diff(x, t); 

% Solve DE with ic
conds = [x(0) == x0, Dx(0) == v0];
sol_with_conds = dsolve(ode, conds);

disp('Solution with Initial Conditions:')
disp(sol_with_conds)

x_t = matlabFunction(sol_with_conds, 'Vars', [t, omega, x0, v0]);

omega_val = 1; 
x0_val = 1;    
v0_val = 0;    

t_vals = linspace(0, 10, 1000);

x_vals = x_t(t_vals, omega_val, x0_val, v0_val);

figure
plot(t_vals, x_vals)
xlabel('Time (t)','interpreter','latex')
ylabel('Displacement (x)','interpreter','latex')
title('Simple Harmonic Oscillator','interpreter','latex')
grid on