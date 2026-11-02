function [x_nom_new, q_ref_new, P_new] = ekf_update_baro_mekf( ...
    x_nom, q_ref, P, z_meas, params)
% EKF_UPDATE_BARO_MEKF  Error-state EKF update for baro altimeter.
%
% Measurement: scalar altitude = z-position.
% Inputs:
%   x_nom (13x1)  - nominal [r(3); v(3); b_g(3); b_a(3); m(1)]
%   q_ref (4x1)   - reference quaternion (Hamilton, scalar-first)
%   P (16x16)     - error-state covariance over [δr; δv; δθ; δb_g; δb_a; δm]
%   z_meas        - baro reading (m)
%   params        - param struct
%
% Outputs: x_nom, q_ref, P after measurement injection (δx reset to zero
% implicitly — never stored).
% L. White, University of Sheffield, May 2026

x_nom = x_nom(:);

%  Observation Jacobian H (1x16) 
% Baro measures r(3), which is error-state index 3.
H = zeros(1, 16);
H(1, 3) = 1;

R_baro = params.baro.noise_std^2;

%  Innovation 
z_corrected = z_meas - params.baro.bias;
z_predicted = x_nom(3);
innov = z_corrected - z_predicted;

%  Kalman gain (16x1) 
S = H * P * H' + R_baro;
K = P * H' / S;


%  Error-state correction δx 
dx = K * innov;     % 16x1

% Unpack the error correction by block
delta_r       = dx(1:3);
delta_v       = dx(4:6);
delta_theta   = dx(7:9);
delta_b_g     = dx(10:12);
delta_b_a     = dx(13:15);
delta_m       = dx(16);

%  Inject error into nominal 
% x_nom layout: [r(3); v(3); b_g(3); b_a(3); m(1)] — indices 1:3, 4:6, 7:9, 10:12, 13
x_nom_new = x_nom;
x_nom_new(1:3)   = x_nom(1:3)   + delta_r;
x_nom_new(4:6)   = x_nom(4:6)   + delta_v;
x_nom_new(7:9)   = x_nom(7:9)   + delta_b_g;
x_nom_new(10:12) = x_nom(10:12) + delta_b_a;
x_nom_new(13)    = x_nom(13)    + delta_m;

% Attitude: multiplicative injection q_ref ← q_ref ⊗ exp_quat(½ δθ)
dq        = exp_quat(delta_theta);
q_ref_new = quat_multiply(q_ref, dq);
q_ref_new = q_ref_new / norm(q_ref_new);

% Mass clamp
x_nom_new(13) = max(x_nom_new(13), params.ms);
x_nom_new(13) = min(x_nom_new(13), params.m0 * 1.10);

% Covariance update (Joseph form) 
I_KH  = eye(16) - K * H;
P_new = I_KH * P * I_KH' + K * R_baro * K';
P_new = (P_new + P_new') / 2;

% Note: the δθ reset Jacobian is approximately I for small corrections
% and is omitted here. This is the standard MEKF simplification and
% is valid as long as |δθ| stays small per update — which it will
% because we apply at 50 Hz to a single scalar measurement.
end


%% LOCAL HELPERS
function dq = exp_quat(delta_theta)
% Build a unit quaternion from a 3-vector rotation in radians.
angle = norm(delta_theta);
if angle > 1e-10
    axis = delta_theta / angle;
    dq   = [cos(angle/2); sin(angle/2) * axis];
else
    dq = [1; 0.5 * delta_theta];
    dq = dq / norm(dq);
end
end

function q_out = quat_multiply(q1, q2)
w1=q1(1); x1=q1(2); y1=q1(3); z1=q1(4);
w2=q2(1); x2=q2(2); y2=q2(3); z2=q2(4);
q_out = [w1*w2 - x1*x2 - y1*y2 - z1*z2;
    w1*x2 + x1*w2 + y1*z2 - z1*y2;
    w1*y2 - x1*z2 + y1*w2 + z1*x2;
    w1*z2 + x1*y2 - y1*x2 + z1*w2];
end