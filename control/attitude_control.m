function [alpha_cmd, beta_cmd, ctrl_out] = attitude_control(x_hat, theta_ref, phi_ref, ...
                                            K_table, ctrl_in, dt, params, omega_body, T_cmd_last)
% ATTITUDE_CONTROL  Gain-scheduled LQI attitude controller.
%
% Runs at 200 Hz. Reads attitude from the MEKF quaternion and mass from the
% 
% MEKF mass state, interpolates pre-computed LQI gains from K_table at the
% current mass, and computes gimbal commands to drive pitch and roll to
% their reference values.
%
% MEKF STATE CONVENTION (Build Plan §8):
%   x_hat = [r(3); v(3); q(4); b_g(3); b_a(3); m]   (17×1)
%   Indices:  1:3    4:6   7:10   11:13   14:16   17
%
%   Attitude quaternion at indices 7:10 (scalar-first: [qw; qx; qy; qz]).
%   Vehicle mass at index 17 (kg).
%
% CONTROL LAW (pitch channel — roll is symmetric):
%
%   e_theta   = theta_hat - theta_ref
%   alpha_cmd = -(K_theta * e_theta + K_q * q_hat + K_i * integral(e_theta))
%
%   e_phi     = phi_hat - phi_ref
%   beta_cmd  = +(K_phi * e_phi + K_p * p_hat + K_i * integral(e_phi))
%
%   Note: beta sign is POSITIVE — the plant gimbal convention for roll
%   (body-x torque, Mx = -T*sin(beta)*ls) requires an inversion relative
%   to pitch to close the loop correctly.
%
% ANGULAR RATE ESTIMATION:
%   Rates (q_hat, p_hat) are estimated by backward-differencing the Euler
%   angle estimates from consecutive MEKF outputs. This avoids derivative
%   kick when the reference steps (we differentiate the measurement, not
%   the error). In the full Simulink system, wire the bias-corrected IMU
%   angular rates here instead for better noise performance.
%
% ANTI-WINDUP:
%   Integrators are clamped to ±i_windup_lim (rad·s) to prevent runaway
%   when the gimbal is saturated (e.g. large reference step or disturbance).
%
% Inputs:
%   x_hat      (17×1)  MEKF state [r; v; q; b_g; b_a; m]
%                      Attitude quaternion at indices 7:10 (scalar-first).
%                      Vehicle mass at index 17 (kg).
%   theta_ref  scalar  commanded pitch angle (rad) from attitude_reference.m
%   phi_ref    scalar  commanded roll  angle (rad) from attitude_reference.m
%   K_table    struct  from lqi_gains() — pre-computed gain lookup table
%   ctrl_in    struct  persistent controller state — init via ctrl_state_init()
%   dt         scalar  timestep (s) — typically 1/200 = 0.005 s
%   params     struct  from plant_params()
%                      Required: params.alpha_max, params.beta_max
%
% Outputs:
%   alpha_cmd  scalar  gimbal pitch command (rad), clamped to ±alpha_max
%   beta_cmd   scalar  gimbal roll  command (rad), clamped to ±beta_max
%   ctrl_out   struct  updated persistent state (pass as ctrl_in next call)
% L. White May 2026



% --- Extract attitude from EKF quaternion (indices 7:10) ---
q_ekf = x_hat(7:10);                             % [qw; qx; qy; qz]
[theta_hat, phi_hat, ~] = quat2euler_local(q_ekf);  % yaw ignored (uncontrolled)

% --- Extract mass from EKF state (index 17) ---
m_est = clamp(x_hat(17), params.ms, params.m0);


if nargin >= 8 && numel(omega_body) == 3
    % Direct body rates — eliminates chattering caused by large angular
    % acceleration at high thrust outrunning the Euler-diff estimator.
    q_rate_hat = omega_body(2);   % body pitch rate (rad/s)
    p_hat      = omega_body(1);   % body roll  rate (rad/s)
else
    q_rate_hat = (theta_hat - ctrl_in.theta_prev) / dt;
    p_hat      = (phi_hat   - ctrl_in.phi_prev)   / dt;
end

% Clamp to reject IMU spikes or initialisation transients
q_rate_hat = clamp(q_rate_hat, -10, 10);
p_hat      = clamp(p_hat,      -10, 10);

% --- Interpolate gains at current mass ---
m_clamped = clamp(m_est, K_table.m_grid(1), K_table.m_grid(end));

K_pitch = interp1(K_table.m_grid, K_table.K_pitch, m_clamped, 'linear');
K_roll  = interp1(K_table.m_grid, K_table.K_roll,  m_clamped, 'linear');

K_theta   = K_pitch(1);   K_q_rate  = K_pitch(2);   K_i_pitch = K_pitch(3);
K_phi     = K_roll(1);    K_p       = K_roll(2);    K_i_roll  = K_roll(3);

% --- Thrust-scheduled scaling ---
% lqi_gains.m synthesises K assuming T = m*g0 (hover thrust).
% Actual control authority b_pitch = T*ls/Ixy scales linearly with T,
% so during coast T is low and the loop is sluggish, during landing
% burn T is high and the loop wants to ring. Rescale gains so the
% closed-loop bandwidth (~sqrt(b*K_theta)) stays roughly constant.
T_design = m_est * params.g0;                  % thrust the gains assume
T_actual = max(T_cmd_last, params.Tmin);       % actual thrust, clamped
K_scale  = min(1.0, T_design / T_actual);   % only reduce gains, never amplify

K_theta   = K_theta   * K_scale;
K_q_rate  = K_q_rate  * sqrt(K_scale);         % damping scales as sqrt(b)
K_i_pitch = K_i_pitch * K_scale;
K_phi     = K_phi     * K_scale;
K_p       = K_p       * sqrt(K_scale);
K_i_roll  = K_i_roll  * K_scale;
% --- Attitude errors ---
e_theta = theta_hat - theta_ref;
e_phi   = phi_hat   - phi_ref;

% --- Compute PD component FIRST for correct anti-windup check ---
alpha_pd = -(K_theta * e_theta + K_q_rate * q_rate_hat);
beta_pd  = (K_phi   * e_phi   + K_p      * p_hat);

% Anti-windup: freeze integrator when PD output alone would saturate
alpha_saturated = abs(alpha_pd) >= params.alpha_max * 0.90;
beta_saturated  = abs(beta_pd)  >= params.beta_max  * 0.90;

i_windup_lim = 0.05;   % rad·s  (was 0.10 — smaller limit reduces integral
                       %          contribution to saturation from 5.2° to 2.3°)

if ~alpha_saturated
    i_theta = clamp(ctrl_in.i_theta + e_theta * dt, -i_windup_lim, i_windup_lim);
else
    i_theta = ctrl_in.i_theta;
end
if ~beta_saturated
    i_phi = clamp(ctrl_in.i_phi + e_phi * dt, -i_windup_lim, i_windup_lim);
else
    i_phi = ctrl_in.i_phi;
end

% --- Control law ---
alpha_cmd = clamp(alpha_pd - K_i_pitch * i_theta, -params.alpha_max, params.alpha_max);
beta_cmd  = clamp(beta_pd  + K_i_roll  * i_phi,   -params.beta_max,  params.beta_max);


% --- Pack persistent state ---
ctrl_out.i_theta      = i_theta;
ctrl_out.i_phi        = i_phi;
ctrl_out.theta_prev   = theta_hat;
ctrl_out.phi_prev     = phi_hat;
ctrl_out.e_theta_prev = e_theta;    % kept for diagnostics
ctrl_out.e_phi_prev   = e_phi;
ctrl_out.m_est        = m_est;      % carried for diagnostics / logging

end



%% LOCAL HELPER

function y = clamp(x, lo, hi)
    y = max(lo, min(hi, x));
end



