classdef studentControllerInterface < matlab.System
    properties (Access = private)
        %% You can add values that you want to store and updae while running your controller.
        % For more information of the supported data type, see
        % https://www.mathworks.com/help/simulink/ug/data-types-supported-by-simulink.html
        t_prev = -1;
        x_hat_prev = [0; 0.00; -pi/3; 0];
        %x_hat_prev = [-0.05; 0.00; 0; 0];
        u_prev = 0;
        theta_d = 0;
        extra_dummy1 = 0;
        extra_dummy2 = 0;
    end
    methods(Access = protected)
        % function setupImpl(obj)
        %    disp("You can use this function for initializaition.");
        % end

        function V_servo = stepImpl(obj, t, p_ball, theta)
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
            y = [p_ball*(20/19.5)-0.0185; theta];
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
            
            theta_d = asin((7 * L / (5 * g * r_arm)) * a_ball_ref);
            u_eq = 0;

            %k_p = 3;
            %theta_d = - k_p * (p_ball - p_ball_ref);
            %theta_saturation = 56 * pi / 180;    
            theta_saturation = 50 * pi / 180;    
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);

            x_ref = [p_ball_ref; v_ball_ref; theta_d; 0];


            % Compute A and B matrices at equilibrium
            A = compute_jacobian_A(x_ref);
            B = compute_jacobian_B();

            % Solve LQR
            Q = diag([800, 0.01, 0.01, 2.5]); % sine -- 0.95, square -- 4
            % Q = diag([1200, 10, 10, 10]);
            R = 0.5;
            % Klqr = lqr(A, B, Q, R);
            

            %Klqr = [10, 25.1525, 13.0233, 2.6315];
            %Klqr = [10, 45, 11, 2.3]; % sine cost: 0.97, square -- 4.4
            %Klqr = [24.5, 33.5, 9.6, 2.7]; % sine cost: 0.85, square --3.5
            %Klqr = [67.08, 66.6, 13.534, 2.7]; % sine cost: , square -- 
            %Klqr = [40, 41.774, 9.123, 1.731]; % sine cost: 0.8081, square --            
            Klqr = [50, 35, 9.123, 1.731]; % sine cost: 0.8081, square -- 
            %Klqr = [40, 45, 9.123, 1.731];
            %Klqr = [5, 20, 0, 0]; % sine cost: 0.8081, square -- 
            

%             if abs(Klqr  * (x_hat - x_ref)) < 0.4
%                 V_servo = u_eq - (Klqr  * (x_hat - x_ref))^(1/3);
%             else
%                 V_servo = u_eq - Klqr  * (x_hat - x_ref);
%             end
            % V_servo = u_eq - Klqr  * (x_hat - x_ref) + 0.4 * sign(Klqr  * (x_hat - x_ref));
            if Klqr  * (x_hat - x_ref) < 0
                V_servo = u_eq - (Klqr  * (x_hat - x_ref)) -0.6;
            else
                V_servo = u_eq - Klqr  * (x_hat - x_ref) +0.6;
            end
            %% nonlinear input control ( works but at a higher energy cost)
            % V_servo = 1*sign(u_eq - Klqr  * (x_hat - x_ref));

            %% saturate V_servo
            lb = -1; % lb = -1 perform better for square
            ub = 1; % ub = 1 perform better for square
            V_servo = min(max(V_servo, lb), ub);


            % % Decide desired servo angle based on simple proportional feedback.
            % k_p = 3;
            % theta_d = - k_p * (p_ball - p_ball_ref);

            % % Make sure that the desired servo angle does not exceed the physical
            % % limit. This part of code is not necessary but highly recommended
            % % because it addresses the actual physical limit of the servo motor.
            % theta_saturation = 56 * pi / 180;    
            % theta_d = min(theta_d, theta_saturation);
            % theta_d = max(theta_d, -theta_saturation);

            % % Simple position control to control servo angle to the desired
            % % position.
            % k_servo = 10;
            % V_servo = k_servo * (theta_d - theta);
            
            % Update class properties if necessary.
            obj.t_prev = t;
            obj.theta_d = theta_d;
            %obj.theta_d = x_hat(3);
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
