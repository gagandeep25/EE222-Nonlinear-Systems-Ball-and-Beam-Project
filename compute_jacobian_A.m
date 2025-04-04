function A = compute_jacobian_A(x_hat)
    p_ball = x_hat(1);
    v_ball = x_hat(2);
    theta  = x_hat(3);
    dtheta = x_hat(4);

    g = 9.81;
    r_arm = 0.0254;
    L = 0.4255;

    a = 5 * g * r_arm / (7 * L);
    b = (5 * L / 14) * (r_arm / L)^2;
    c = (5 / 7) * (r_arm / L)^2;

    A = [ 0, 1, 0, 0;
          0, 0, a * cos(theta) + 2 * b * dtheta^2 * cos(theta) * sin(theta) - 2 * c * p_ball * cos(theta) * sin(theta), -2 * b * dtheta * cos(theta)^2 + 2 * c * p_ball * dtheta * cos(theta)^2;
          0, 0, 0, 1;
          0, 0, 0, -1/0.025];
end

