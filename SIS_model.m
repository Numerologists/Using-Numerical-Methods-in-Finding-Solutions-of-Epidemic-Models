% Parameters
r = 0.1;     % transmission rate
gamma = 4.5;    % recovery rate
S0 = 25;       % initial susceptible population
I0 = 10;         % initial infectious population
h = 0.1;        % time step
num_steps = 50;

% Initialize arrays to store results
t_values = zeros(num_steps + 1, 1);
s_values = zeros(num_steps + 1, 1);
i_values = zeros(num_steps + 1, 1);

% Initial conditions
s_values(1) = S0;
i_values(1) = I0;

% Perform RK4 integration
for n = 1:num_steps
    t_values(n+1) = t_values(n) + h;
    y = [s_values(n); i_values(n)];
    y_next = rk4_step(@sis_model, t_values(n), y, h, r, gamma);
    s_values(n+1) = y_next(1);
    i_values(n+1) = y_next(2);
end

% Calculate total population K = S + I
K_values = s_values + i_values;

% Plot I vs K with time
figure;
%plot(t_values, K_values, 'b-', 'LineWidth', 2, 'DisplayName', 'Total Population (K)');
%hold on;
plot(t_values, i_values, 'r-', 'LineWidth', 2, 'DisplayName', 'Infected Population (I)');
hold on;
plot(t_values, s_values, 's-', 'LineWidth', 2, 'DisplayName', 'Susceptible Population  (S)');
xlabel('Time');
ylabel('Population');
title(' when rK is less than gamma');
legend('show');
grid on;

% Function to define the SIS model equations
function dydt = sis_model(t, y, r, gamma)
    s = y(1);
    i = y(2);
    dsdt = -r * s * i + gamma * i;
    didt = r * s * i - gamma * i;
    dydt = [dsdt; didt];
end

% Function to perform RK4 integration
function y_next = rk4_step(func, t, y, h, r, gamma)
    k1 = func(t, y, r, gamma);
    k2 = func(t + 0.5 * h, y + 0.5 * h * k1, r, gamma);
    k3 = func(t + 0.5 * h, y + 0.5 * h * k2, r, gamma);
    k4 = func(t + h, y + h * k3, r, gamma);
    y_next = y + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
end
