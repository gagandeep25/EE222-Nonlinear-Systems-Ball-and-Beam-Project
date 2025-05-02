classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        x_hat_prev = [0; 0.00; -pi/3; 0]; % initial condition for hardware
        %x_hat_prev = [-0.05; 0.00; 0; 0];
        u_prev = 0;
        theta_d = 0;
        extra_dummy1 = 0;
        extra_dummy2 = 0;
    end
    methods(Access = protected)
        function V_servo = stepImpl(obj, t, p_ball, theta)
           
        % LQR controller implementation
        % This is the main function called every iteration. You have to implement
        % the controller in this function, bu you are not allowed to
        % change the signature of this function. 
        % Input arguments:
        %   t: current time
        %   p_ball: position of the ball provided by the ball position sensor (m)
        %
        %   theta: servo motor angle provided by the encoder of the motor (rad)
        % Output:
        %   V_servo: voltage to the servo input.        
            %% Sample Controller: Simple Proportional Controller
            t_prev = obj.t_prev;
            theta_d = obj.theta_d;
            x_hat_prev = obj.x_hat_prev;
            y = [p_ball*(20/19.5)-0.0185; theta]; % adjustment to the sensor output to account for measurement errors
            %y = [p_ball; theta];
            u_prev = 0;

             % System parameters
             g = 9.81;
             r_arm = 0.0254;
             L = 0.4255;
             K = 1.5;
                      
            % Extract reference trajectory at the current timestep.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
            x_op = [p_ball_ref, v_ball_ref, 0, 0];

            % state_estimate -- luenberger observer
            x_hat = luenberger_observer(t-t_prev, x_hat_prev, y, u_prev, x_op);
            
            theta_d = asin((7 * L / (5 * g * r_arm)) * a_ball_ref); % theta_d estimation from reference trajectory
            u_eq = 0;

            theta_saturation = 50 * pi / 180;   % limiting the angle to +/- 50 degrees 
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            x_ref = [p_ball_ref; v_ball_ref; theta_d; 0];


            % % Compute A and B matrices at equilibrium
            % Time varying implementation
            % A = compute_jacobian_A(x_ref);
            % B = compute_jacobian_B();

            % % Solve LQR
            % Q = diag([800, 0.01, 0.01, 2.5]); % sine -- 0.95, square -- 4
            % Q = diag([1200, 10, 10, 10]);
            % R = 0.5;
            % Klqr = lqr(A, B, Q, R);
            
            % Time invariant Klqr calculated at x_ref = [0, 0, 0, 0]
            %Klqr = [40, 41.774, 9.123, 1.731]; % used in simulation         
            Klqr = [50, 35, 9.123, 1.731]; % used for hardware
            
            if Klqr  * (x_hat - x_ref) < 0
                V_servo = u_eq - (Klqr  * (x_hat - x_ref)) -0.6;
            else
                V_servo = u_eq - Klqr  * (x_hat - x_ref) +0.6;
            end

            % % experimentally reduced the total magnitude of V_servo; only
            % for hardware
            V_servo = V_servo * 0.8;
            % % saturate V_servo to +/-1 V
            lb = -1;
            ub = 1; 
            V_servo = min(max(V_servo, lb), ub);
            
            % Update class properties if necessary.
            obj.t_prev = t;
            obj.theta_d = theta_d;
            obj.x_hat_prev = x_hat;
            obj.u_prev = V_servo;
        end
    end
    
    methods(Access = public)
        % Used this for matlab simulation script. fill free to modify it as
        % however you want.
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)        
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
    
end
