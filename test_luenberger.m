clc; clear; close all;

% Simulation Parameters
delta_t = 0.001;  % Time step (1 ms)
T = 5;            % Total simulation time (seconds)
time = 0:delta_t:T;

% Initial Conditions
x_true = [0.05; 0; 0.1; 0];   % [p_ball, v_ball, theta, dtheta]
x_hat = [0; 0; 0.1; 0];         % Initial observer estimate
u_func = @(t) 0.1 * (t > 1);  % Control input: step at t=1s

% Data Storage
x_true_log = zeros(4, length(time));
x_hat_log = zeros(4, length(time));

% Run Simulation
for i = 1:length(time)
    t = time(i);
    u = u_func(t);
    
    % True system update (Euler method for simplicity)
    dx_true = ball_and_beam_dynamics(x_true, u);
    x_true = x_true + delta_t * dx_true;
    
    % Observer update
    y_measured = x_true([1, 3]);  % Only measure p_ball and theta
    x_hat = luenberger_observer(delta_t, x_hat, y_measured, u);
    
    % Store results
    x_true_log(:, i) = x_true;
    x_hat_log(:, i) = x_hat;
end

% Plot results
figure;
titles = {'Ball Position (m)', 'Ball Velocity (m/s)', 'Beam Angle (rad)', 'Beam Angular Velocity (rad/s)'};
for i = 1:4
    subplot(2,2,i);
    plot(time, x_true_log(i,:), 'b', 'LineWidth', 1.5); hold on;
    plot(time, x_hat_log(i,:), 'r--', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel(titles{i});
    legend('True', 'Estimated');
    grid on;
end
sgtitle('Luenberger Observer Performance');