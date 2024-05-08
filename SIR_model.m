

% Parameters
r = 0.0178;     % transmission rate (replaced beta with r)
alpha = 2.73; % recovery rate
S0 = 254;    % initial susceptible population
I0 = 7;    % initial infectious population
h = 0.1;     % time step
num_steps = 60;

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
    y_next = rk4_step(@sir_model, t_values(n), y, h, r, alpha);
    s_values(n+1) = y_next(1);
    i_values(n+1) = y_next(2);
end

% Plot the results
figure;
plot(t_values, s_values, 'b-', 'LineWidth', 2, 'DisplayName', 'Susceptible');
hold on;
%plot(t_values, i_values, 'r-', 'LineWidth', 2, 'DisplayName', 'Infectious');
xlabel('Time');
ylabel('Population');
title('');
legend('show');
grid on;
% Function to define the SIR model equations
function dydt = sir_model(t, y, r, alpha)
    s = y(1);
    i = y(2);
    dsdt = -r * s * i;
    didt = r * s * i - alpha * i;
    dydt = [dsdt; didt];
end

% Function to perform RK4 integration
function y_next = rk4_step(func, t, y, h, r, alpha)
    k1 = func(t, y, r, alpha);
    k2 = func(t + 0.5 * h, y + 0.5 * h * k1, r, alpha);
    k3 = func(t + 0.5 * h, y + 0.5 * h * k2, r, alpha);
    k4 = func(t + h, y + h * k3, r, alpha);
    y_next = y + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
end
