function [T_cmd, n_hat_cmd, tracking_error] = guidance_track_reference(t, x_hat, reference, params)
% GUIDANCE_TRACK_REFERENCE  Mode A reference-tracking guidance.
%
% Looks up the pre-computed ICLOCS reference trajectory, extracts a
% feedforward thrust direction and magnitude at the current time, and
% adds a PD correction proportional to position and velocity error.
% Used proactively during ascent and as a fallback if Mode B rejects
% three consecutive solutions.
%
% INPUTS:
%   t         - current mission time                        [s]       scalar
%   x_hat     - EKF state estimate                         [mixed]   17x1
%                 x_hat(1:3)  = inertial position r        [m]
%                 x_hat(4:6)  = inertial velocity v        [m/s]
%                 x_hat(7:10) = attitude quaternion q       [-]
%                 x_hat(11:13)= gyro  bias b_g             [rad/s]
%                 x_hat(14:16)= accel bias b_a             [m/s²]
%                 x_hat(17)   = vehicle mass m             [kg]
%   reference - struct loaded from reference_trajectory.mat
%                 .t     (4000x1) time stamps              [s]
%                 .x     (4000x13) state history           [mixed]
%                 .u     (4000x3) control history          [N, deg, deg]
%                 .phase (4000x1) phase index              [1–4]
%   params    - vehicle params struct (plant_params.m)
%
% OUTPUTS:
%   T_cmd          - commanded thrust magnitude            [N]       scalar
%   n_hat_cmd      - commanded thrust unit direction       [-]       3x1 inertial
%   tracking_error - struct with fields:
%                      .pos_err    position error norm     [m]       scalar
%                      .vel_err    velocity error norm     [m/s]     scalar
%                      .r_err_vec  position error vector   [m]       3x1
%                      .v_err_vec  velocity error vector   [m/s]     3x1
%                      .phase_id   current phase index     [-]       scalar
%
% CALLS:     alpha_beta_to_nhat.m  (local helper, defined below)
% CALLED BY: guidance subsystem block in vtvl_sim.slx
% L. White, University of Sheffield, May 2026
%%  STEP 1 — Interpolate reference trajectory at current time t
% Pre-processing: remove duplicate time entries at phase boundaries.
% ICLOCS stores the last node of one phase and the first node of the
% next at identical timestamps. interp1 requires strictly unique
% sample points, so duplicates are removed here (keeping the first
% occurrence via 'stable', the phase-end value).
[t_clean, ia]   = unique(reference.t,    'stable');
x_clean         = reference.x(ia, :);
u_clean         = reference.u(ia, :);
phase_clean     = reference.phase(ia, :);
% t_clean is now 3997 points (4000 minus 3 duplicate boundaries)

r_ref = interp1(t_clean, x_clean(:, 1:3), t, 'linear', 'extrap')';
% r_ref: inertial reference position [m], 3x1

v_ref = interp1(t_clean, x_clean(:, 4:6), t, 'linear', 'extrap')';
% v_ref: inertial reference velocity [m/s], 3x1

T_ref = interp1(t_clean, u_clean(:, 1), t, 'linear', 'extrap');
% T_ref: reference thrust magnitude [N], scalar

q_ref = interp1(t_clean, x_clean(:, 8:11), t, 'linear', 'extrap')';
% q_ref: reference attitude quaternion [-]  4x1

alpha_ref_deg = interp1(t_clean, u_clean(:, 2), t, 'linear', 'extrap');
beta_ref_deg  = interp1(t_clean, u_clean(:, 3), t, 'linear', 'extrap');
% alpha_ref_deg, beta_ref_deg: gimbal angles [deg] — converted in Step 5


%%  STEP 2- Extract current state from MEKF estimate
% Pull position, velocity, and mass out of the 17-element
%       x_hat vector by index.
r_hat = x_hat(1:3);    % inertial position estimate  [m]     3x1
v_hat = x_hat(4:6);    % inertial velocity estimate  [m/s]   3x1
m_hat = x_hat(17);     % mass estimate               [kg]    scalar

% Guard: clamp mass to physical bounds to protect Step 5 division.
m_hat = max(m_hat, params.ms);

%%  STEP 3 - Compute tracking errors
% Subtract the MEKF estimate from the reference to get the position and velocity errors the PD correction must null.
e_r = r_ref - r_hat;    % position error vector   [m]     3x1
e_v = v_ref - v_hat;    % velocity error vector   [m/s]   3x1

%%  STEP 4 - PD correction and thrust command
% Combine the feedforward thrust direction from the reference with a PD 
% correction that nulls the tracking errors from Step 4. Recover a thrust
%  magnitude and unit direction.
%
% Get feedforward thrust direction from reference
% Convert reference gimbal angles (degrees) to a unit thrust
% direction vector in the inertial frame using the helper below.
% This is the direction the ICLOCS trajectory intended to point.
n_ref = alpha_beta_to_nhat(alpha_ref_deg, beta_ref_deg, q_ref);
% Uses reference attitude — feedforward is pure ICLOCS intent, 
% not contaminated by current attitude error
% n_ref: reference thrust unit direction [inertial frame]  3x1

% Feedforward acceleration
% Reconstruct the acceleration the reference trajectory intended
% to produce at this moment. Dividing by m_hat (not reference mass)
% means if the vehicle is heavier or lighter than expected, the
% feedforward term automatically adjusts.
a_ff = (T_ref / m_hat) * n_ref;    % feedforward acceleration [m/s²]  3x1

% PD correction acceleration
% The proportional term (Kp * e_r) corrects slow position drift.
% The derivative term (Kd * e_v) damps oscillations before they grow.
% Both produce accelerations [m/s²] — gains have units s^-2 and s^-1.
a_corr = params.modeA.Kp * e_r + params.modeA.Kd * e_v; %3x1  [m/s²]

% Total commanded acceleration
a_cmd = a_ff + a_corr;             % total acceleration command [m/s²]  3x1

% Recover thrust magnitude and direction
% T_cmd = m_hat * ||a_cmd||, then clamped to physical actuator limits.
% If a_cmd is somehow zero (e.g. t=0 with perfect tracking), the
% fallback points thrust straight up — vehicle holds altitude rather
% than going to zero thrust and falling.
T_cmd = max(params.Tmin, min(params.Tmax, m_hat * norm(a_cmd)));
% T_cmd: scalar thrust magnitude [N], clamped to [Tmin, Tmax]

if norm(a_cmd) > 1e-6
    n_hat_cmd = a_cmd / norm(a_cmd);    % unit thrust direction  3x1
else
    n_hat_cmd = [0; 0; 1];             % fallback: point straight up
end
%%  STEP 6 — Populate tracking_error struct for logging
% Pack all diagnostic quantities into a struct that the
%  Simulink logger and verification scripts can read.
tracking_error.pos_err   = norm(e_r);   % scalar position error norm  [m]
tracking_error.vel_err   = norm(e_v);   % scalar velocity error norm  [m/s]
tracking_error.r_err_vec = e_r;         % 3x1 position error vector   [m]
tracking_error.v_err_vec = e_v;         % 3x1 velocity error vector   [m/s]
tracking_error.phase_id  = interp1(t_clean, double(phase_clean), ...
    t, 'nearest', 'extrap');
% phase_id: integer phase label (1-4) at current time
end

%%  LOCAL HELPER — alpha_beta_to_nhat
function n_hat = alpha_beta_to_nhat(alpha_deg, beta_deg, q)
% ALPHA_BETA_TO_NHAT  Convert gimbal angles to inertial thrust direction.
%
% Converts reference gimbal angles (stored in degrees in reference.u)
% into a unit thrust direction vector in the inertial frame, using the
% current attitude quaternion to rotate from body to inertial frame.
%
% INPUTS:
%   alpha_deg - gimbal pitch angle  [deg]   scalar
%   beta_deg  - gimbal yaw angle    [deg]   scalar
%   q         - attitude quaternion [-]     4x1  [q0;q1;q2;q3]
%
% OUTPUTS:
%   n_hat     - thrust unit direction [inertial frame]  3x1

% Convert to radians
alpha = alpha_deg * pi/180;
beta  = beta_deg  * pi/180;

% Thrust direction in body frame (matches plant_dynamics.m Step 4)
n_body = [ sin(alpha);          
          -sin(beta);           
           cos(alpha)*cos(beta)];
% Normalise body vector (should already be unit but guard against
% numerical error near alpha=beta=0)
n_body = n_body / norm(n_body);

% Rotation matrix body -> inertial from quaternion
q0=q(1); q1=q(2); q2=q(3); q3=q(4);
R_bi = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
    2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
    2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];

% Rotate to inertial frame
n_hat = R_bi * n_body;
n_hat = n_hat / norm(n_hat);    % renormalise after rotation
end

