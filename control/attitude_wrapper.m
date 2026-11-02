function [alpha_cmd, beta_cmd] = attitude_wrapper(x_hat, n_hat_cmd, gyro_meas, t, T_cmd_last)
% ATTITUDE_WRAPPER  Converts a thrust direction command into attitude angle
% commands (alpha, beta) via the attitude controller.
%
% Inputs:
%   x_hat       - full state estimate vector (17x1)
%   n_hat_cmd   - desired thrust direction unit vector (3x1)
%   gyro_meas   - raw gyroscope measurement (3x1) [rad/s]
%   t           - current time [s], used to detect simulation restart
%   T_cmd_last  - thrust magnitude from previous step [N]
%
% Outputs:
%   alpha_cmd   - commanded pitch angle [rad]
%   beta_cmd    - commanded roll angle [rad]
% L. White May 2026

persistent ctrl_state params_cached

% Load params from base workspace once and cache for subsequent calls
if isempty(params_cached)
    params_cached = evalin('base', 'params');
end
params = params_cached;

% Initialise controller state at startup or if time has reset
if isempty(ctrl_state) || t < 0.005
    ctrl_state = ctrl_state_init(params.m0);
end

dt = 0.005;   % fixed attitude control timestep [s]

% Ensure all inputs are column vectors
x_hat_col = x_hat(:);
n_hat_col = n_hat_cmd(:);
gyro_col  = gyro_meas(:);

% Convert thrust direction to reference pitch and roll angles
[theta_ref, phi_ref] = attitude_reference(n_hat_col);

% Subtract estimated gyro bias to get true body-frame angular rate
b_g        = x_hat_col(11:13);
omega_body = gyro_col - b_g;

% Run attitude controller and update persistent state
[alpha_cmd, beta_cmd, ctrl_state] = attitude_control(x_hat, theta_ref, phi_ref, ...
   params.K_table, ctrl_state, dt, params, omega_body, T_cmd_last);
end