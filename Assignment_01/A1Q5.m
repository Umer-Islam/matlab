% Define and plot a function using the symbolic math toolbox.
syms x
f = exp(x)-3;

%  derivative
df = diff(f, x);

x_values = linspace(-3, 3, 100);
f_values = double(subs(f, x, x_values));
df_values = double(subs(df, x, x_values));

% Plot 
figure;
hold on;
plot(x_values, f_values, 'LineWidth', 2, 'DisplayName', 'f(x)', 'Color', '#FFA500');

plot(x_values, df_values, 'LineWidth', 2, 'DisplayName', "f'(x)", 'Color', '#00FF00'); 
legend('Location', 'Best');
title('Function and Its Derivative');
xlabel('x');
ylabel('y');
grid on;
hold off;
