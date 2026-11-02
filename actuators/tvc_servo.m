function [alpha_act, beta_act, act_out] = tvc_servo(alpha_cmd, beta_cmd, act_in, dt, params)
% TVC_SERVO  Second-order servo model for the TVC gimbal actuators.
%
% Models each gimbal axis as a second-order system with:
%   - Natural frequency  wn   = 25 rad/s  (bandwidth ~4 Hz)
%   - Damping ratio      zeta = 0.7       (slightly underdamped — fast but stable)
%   - Rate limit         60 deg/s         (slew-rate saturation)
%   - Position limit     ±alpha_max / ±beta_max  (from params)
%
% NOTE:
%   The rate limiter tracks the PREVIOUS RATE-LIMITED COMMAND (alpha_cmd_rl),
%   NOT the servo output.  Tracking servo output causes a silent failure:
%   the limiter only advances one lag-epsilon per step (output ≈ cmd_rl - ε),
%   reducing effective slew to ~5 deg/s instead of 60 deg/s.
%
% Discretised via Tustin (bilinear) transform at the controller rate.
%
% Inputs:
%   alpha_cmd  (rad)    commanded gimbal pitch angle from attitude_control()
%   beta_cmd   (rad)    commanded gimbal roll  angle from attitude_control()
%   act_in     struct   persistent actuator state from act_state_init()
%   dt         (s)      timestep (~0.005 s at 200 Hz)
%   params     struct   from plant_params() — uses alpha_max, beta_max
%
% Outputs:
%   alpha_act  (rad)    actual gimbal pitch angle delivered to plant
%   beta_act   (rad)    actual gimbal roll  angle delivered to plant
%   act_out    struct   updated persistent state
% L. White, University of Sheffield, May 2026

% --- Servo parameters ---
wn       = 25.0;         % natural frequency (rad/s)
zeta     = 0.70;         % damping ratio
rate_lim = 60*pi/180;    % max slew rate (rad/s) = 60 deg/s

% --- Continuous-time state-space (per axis) ---
Ac = [0,       1;
     -wn^2,   -2*zeta*wn];
Bc = [0; wn^2];
Cc = [1, 0];

% --- Tustin discretisation ---
I2  = eye(2);
tmp = I2 - (dt/2)*Ac;
Ad  = tmp \ (I2 + (dt/2)*Ac);
Bd  = tmp \ (dt * Bc);

% --- Pitch axis ---
[alpha_act, s_alpha, alpha_cmd_rl] = servo_axis( ...
    alpha_cmd, act_in.s_alpha, act_in.alpha_cmd_rl, ...
    rate_lim, params.alpha_max, Ad, Bd, Cc, dt);

% --- Roll axis ---
[beta_act, s_beta, beta_cmd_rl] = servo_axis( ...
    beta_cmd, act_in.s_beta, act_in.beta_cmd_rl, ...
    rate_lim, params.beta_max, Ad, Bd, Cc, dt);

% --- Pack state ---
act_out.s_alpha      = s_alpha;
act_out.s_beta       = s_beta;
act_out.alpha_cmd_rl = alpha_cmd_rl;  % rate-limited command trajectory
act_out.beta_cmd_rl  = beta_cmd_rl;
act_out.alpha_prev   = alpha_act;     % actual servo output (legacy field)
act_out.beta_prev    = beta_act;

end


function [y, s_new, cmd_rl] = servo_axis(cmd, s, cmd_rl_prev, rate_lim, pos_lim, Ad, Bd, Cc, dt)
% One gimbal axis: rate-limit on cmd_rl trajectory, then 2nd-order dynamics.

    % Rate-limit relative to previous rate-limited command (NOT servo output)
    delta     = cmd - cmd_rl_prev;
    max_delta = rate_lim * dt;
    cmd_rl    = cmd_rl_prev + clamp(delta, -max_delta, max_delta);

    % Second-order servo dynamics
    s_new = Ad * s + Bd * cmd_rl;

    % Position output + hard stop
    y        = clamp(Cc * s_new, -pos_lim, pos_lim);
    s_new(1) = y;   % write back to keep state consistent with saturation
end


function y = clamp(x, lo, hi)
    y = max(lo, min(hi, x));
end