function [x_nom_new, q_ref_new, P_new] = ekf_predict_mekf( ...
        x_nom, q_ref, P, accel_meas, gyro_meas, dt, params, T_cmd)
% EKF_PREDICT_MEKF  Error-state EKF predict step.
%
% Propagates the nominal state forward and updates the 16×16 error-state
% covariance using the correct linearised Jacobians for body-frame δθ.
% L. White, University of Sheffield, May 2026

x_nom      = x_nom(:);
accel_meas = accel_meas(:);
gyro_meas  = gyro_meas(:);

%  Unpack nominal 
r   = x_nom(1:3);
v   = x_nom(4:6);
b_g = x_nom(7:9);
b_a = x_nom(10:12);
m   = x_nom(13);

%  Bias-corrected IMU 
omega = gyro_meas  - b_g;       % estimated angular rate, body frame
f     = accel_meas - b_a;       % estimated specific force, body frame

%  Rotation from current q_ref 
R = quat2rotm_local(q_ref);     % body → inertial

%  Strapdown propagation of r, v 
g_vec      = [0; 0; -params.g0];
f_inertial = R * f;
a          = f_inertial + g_vec;
r_new = r + v*dt + 0.5*a*dt^2;
v_new = v + a*dt;

%  Multiplicative propagation of q_ref 
% q_ref(k+1) = q_ref(k) ⊗ exp_quat(½ ω dt)
% This keeps q_ref exactly on the unit sphere by construction.
angle = norm(omega) * dt;
if angle > 1e-10
    axis = omega / norm(omega);
    dq   = [cos(angle/2); sin(angle/2) * axis];
else
    dq = [1; 0.5 * omega * dt];
    dq = dq / norm(dq);
end
q_ref_new = quat_multiply(q_ref, dq);
q_ref_new = q_ref_new / norm(q_ref_new);

%  Biases and mass 
b_g_new = b_g;                  % random-walk: mean unchanged in predict
b_a_new = b_a;
T_mag   = max(T_cmd, 0);
m_new   = m - (T_mag / params.c) * dt;
m_new   = max(m_new, params.ms);

%  Assemble new nominal 
x_nom_new = [r_new; v_new; b_g_new; b_a_new; m_new];

%% ERROR-STATE F MATRIX (16×16) — the key MEKF math
% State order:  [δr(3); δv(3); δθ(3); δb_g(3); δb_a(3); δm(1)]
%
% Linearised continuous-time dynamics:
%   δṙ = δv
%   δv̇ = -[R f]× δθ - R δb_a + noise
%   δθ̇ = -[ω]× δθ - δb_g + noise
%   δb_g, δb_a, δm = random walks
%
% In discrete time, with first-order Euler:
F = eye(16);

% Position from velocity error
F(1:3, 4:6) = eye(3) * dt;

% Velocity from attitude error
F(4:6, 7:9) = -R * skew(f) * dt;

% Velocity from accel-bias error
F(4:6, 13:15) = -R * dt;

% Attitude from attitude (rotation drift)
F(7:9, 7:9) = eye(3) - skew(omega) * dt;

% Attitude from gyro-bias error
F(7:9, 10:12) = -eye(3) * dt;

% Mass row uncoupled; F(16,16) = 1 already from eye(16).

%% PROCESS NOISE Q (16×16)
p = params.ekf;
q_accel = p.accel_noise_std.^2;
q_gyro  = p.gyro_noise_std.^2;
q_m     = p.mass_process_noise_std^2;

Q = zeros(16, 16);
Q(4:6,   4:6)   = diag(q_accel) * dt;   % accel noise drives δv
Q(7:9,   7:9)   = diag(q_gyro)  * dt;   % gyro noise drives δθ
Q(10:12, 10:12) = diag(p.gyro_bias_rw.^2)  * dt;
Q(13:15, 13:15) = diag(p.accel_bias_rw.^2) * dt;
Q(16, 16)       = q_m * dt;

% DIAGNOSTIC: freeze biases (zero out the bias rows only) 
Q(10:12, 10:12) = zeros(3,3);
Q(13:15, 13:15) = zeros(3,3);

% Covariance propagation 
P_new = F * P * F' + Q;
P_new = (P_new + P_new') / 2;
end

%% LOCAL HELPERS

function S = skew(v)
% SKEW  Skew-symmetric matrix of a 3-vector, so that skew(a)*b = a × b.
S = [   0,  -v(3),  v(2);
      v(3),    0,  -v(1);
     -v(2),  v(1),    0];
end

function R = quat2rotm_local(q)
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [1-2*(qy^2+qz^2),   2*(qx*qy-qw*qz),   2*(qx*qz+qw*qy);
     2*(qx*qy+qw*qz),   1-2*(qx^2+qz^2),   2*(qy*qz-qw*qx);
     2*(qx*qz-qw*qy),   2*(qy*qz+qw*qx),   1-2*(qx^2+qy^2)];
end

function q_out = quat_multiply(q1, q2)
w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
q_out = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
         w1*x2 + x1*w2 + y1*z2 - z1*y2;
         w1*y2 - x1*z2 + y1*w2 + z1*x2;
         w1*z2 + x1*y2 - y1*x2 + z1*w2];
end