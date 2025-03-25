function x_hat = luenberger_observer(delta_t, x_hat_prev, y, u)
    % x_hat_prev -- 4x1 vector
    % y -- 2x1 vector

    LO = zeros(4,2); %Luenberger gain
    LO = [[  25.1320,  -0.5331]; [ 151.9416,  -6.4568]; [  -1.2948, -10.1320]; [  32.6879, 619.3871]];
    C = [[1, 0, 0, 0]; [0, 0, 1, 0]]; %selects position and theta

    g = 9.81;
    r_arm = 0.0254;
    L = 0.4255;

    a = 5 * g * r_arm / (7 * L);
    b = (5 * L / 14) * (r_arm / L)^2;
    c = (5 / 7) * (r_arm / L)^2;

    dx_hat = zeros(4, 1);
    x_hat = zeros(4, 1);

    % x_hat dynamics
    dx_hat(1) = x_hat_prev(2);
    dx_hat(2) = a * sin(x_hat_prev(3)) - b * x_hat_prev(4)^2 * cos(x_hat_prev(3))^2 + c * x_hat_prev(1) * x_hat_prev(4)^2 * cos(x_hat_prev(3))^2;
    dx_hat(3) = x_hat_prev(4);
    K = 1.5;
    tau = 0.025;
    dx_hat(4) = (- x_hat_prev(4) + K * u) / tau; 

    % Luenberger part
    dx_hat = dx_hat + LO * (y - C*x_hat_prev);

    % x_hat update step
    x_hat = x_hat_prev + dx_hat * delta_t;
end
