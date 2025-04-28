% PID CONTROLLER
classdef studentControllerInterface < matlab.System
    properties (Access = private)
        % Existing properties
        t_prev = -1;
        x_hat_prev = [0; 0.00; -pi/3; 0];
        theta_d = 0;
        extra_dummy1 = 0;
        extra_dummy2 = 0;
        % New properties for PID control
        e_prev = 0;
        e_integral = 0;
        u_prev = 0;
    end
    methods(Access = protected)
        function V_servo = stepImpl(obj, t, p_ball, theta)
            % Get the reference trajectory at time t.
            [p_ball_ref, v_ball_ref, a_ball_ref] = get_ref_traj(t);
             %y = [p_ball; theta];
            y = [p_ball*(20/19.5)-0.0185; theta];

            x_op = [p_ball_ref, v_ball_ref, 0, 0];

            % state_estimate -- luenberger observer
            x_hat = luenberger_observer(t-obj.t_prev, obj.x_hat_prev, y, obj.u_prev, x_op);
            error = (x_hat(1) - p_ball_ref) + 0.5*(x_hat(2) - v_ball_ref); 
            obj.x_hat_prev = x_hat;
            
            % Compute the error between the actual ball position and its reference.
            % error = p_ball - p_ball_ref;
            
            %% PID Gains
            % Adjust these gains to tune the controller behavior:
            % k_p: Proportional gain (increases response speed but too high may lead to instability). 15
            % k_i: Integral gain (reduces steady-state error; too high may cause overshoot). 0.75
            % k_d: Derivative gain (dampens the response; too high may lead to noise sensitivity). 40
            k_p = 2;   % Proportional gain (ADJUST) 80, 60
            k_i = 0.4;   % Integral gain (ADJUST) 
            k_d = 0; % Derivative gain (ADJUST)40
            
            % Compute the time step (dt) and handle the first iteration.
            if obj.t_prev < 0
                dt = 0;
                derivative = 0;
                obj.e_integral = 0;
            else
                dt = t - obj.t_prev;
                % Avoid division by zero if dt is zero or negative.
                if dt <= 0
                    dt = eps;
                end
                % Compute the derivative (rate of change of error).
                derivative = (error - obj.e_prev) / dt;
                % Accumulate the error over time for the integral term.
                obj.e_integral = obj.e_integral + error * dt;
            end
            
            % Compute the desired servo angle using the PID control law.
            % A negative sign is used to match the original control direction.
            theta_d = - (k_p * error + k_i * obj.e_integral + k_d * derivative);
            
            %% Saturation Limits
            % Make sure that the desired servo angle does not exceed the physical
            % limit. This part of code is not necessary but highly recommended
            % because it addresses the actual physical limit of the servo motor.
            theta_saturation = 50 * pi / 180;    
            theta_d = min(theta_d, theta_saturation);
            theta_d = max(theta_d, -theta_saturation);
            
            %% Servo Voltage Gain
            % Simple position control to control servo angle to the desired
            % position.
            % Adjust the servo gain (k_servo) to scale the voltage command appropriately.12
            k_servo = 3;  % Servo voltage gain (ADJUST) 15
            V_servo = k_servo * (theta_d - theta);
            lb = -1;
            ub = 1;
            if V_servo < 0
                V_servo = V_servo - 0.6;
            else
                V_servo = V_servo + 0.6;
            end
            V_servo = min(max(V_servo, lb), ub);
            % Update properties for the next iteration.
            obj.t_prev = t;
            obj.theta_d = theta_d;
            obj.e_prev = error;
            obj.u_prev = V_servo;
        end
    end
    
    methods(Access = public)
        function [V_servo, theta_d] = stepController(obj, t, p_ball, theta)
            V_servo = stepImpl(obj, t, p_ball, theta);
            theta_d = obj.theta_d;
        end
    end
end
