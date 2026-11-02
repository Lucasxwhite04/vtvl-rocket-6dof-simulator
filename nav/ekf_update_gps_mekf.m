function [x_nom_new, q_ref_new, P_new] = ekf_update_gps_mekf( ...
    x_nom, q_ref, P, r_meas, v_meas, params)
% EKF_UPDATE_GPS_MEKF  Error-state EKF update for GPS pos + vel.
%
% Measurement: [r(3); v(3)] = 6-D vector.
% Updates 16-state error: [δr; δv; δθ; δb_g; δb_a; δm]
% Then injects the correction into x_nom and q_ref.
% L. White, University of Sheffield, May 2026

x_nom  = x_nom(:);
r_meas = r_meas(:);
v_meas = v_meas(:);

%  Observation matrix H (6x16) 
% Measures δr (rows 1-3) and δv (rows 4-6) of the error state.
H = zeros(6, 16);
H(1:3, 1:3) = eye(3);    % δr
H(4:6, 4:6) = eye(3);    % δv

R_gps = diag([params.gps.pos_noise_std.^2; params.gps.vel_noise_std.^2]);

%  Innovation 
z       = [r_meas; v_meas];
z_pred  = [x_nom(1:3); x_nom(4:6)];
innov   = z - z_pred;

%  Kalman gain (16x6) 
S = H * P * H' + R_gps;
K = P * H' / S;

%  Error-state correction δx 
dx = K * innov;     % 16x1

delta_r       = dx(1:3);
delta_v       = dx(4:6);
delta_theta   = dx(7:9);
delta_b_g     = dx(10:12);
delta_b_a     = dx(13:15);
delta_m       = dx(16);

%  Inject into nominal 
x_nom_new = x_nom;
x_nom_new(1:3)   = x_nom(1:3)   + delta_r;
x_nom_new(4:6)   = x_nom(4:6)   + delta_v;
x_nom_new(7:9)   = x_nom(7:9)   + delta_b_g;
x_nom_new(10:12) = x_nom(10:12) + delta_b_a;
x_nom_new(13)    = x_nom(13)    + delta_m;

% Attitude: multiplicative
dq        = exp_quat(delta_theta);
q_ref_new = quat_multiply(q_ref, dq);
q_ref_new = q_ref_new / norm(q_ref_new);

% Mass clamp
x_nom_new(13) = max(x_nom_new(13), params.ms);
x_nom_new(13) = min(x_nom_new(13), params.m0 * 1.10);

%  Covariance update (Joseph form) 
I_KH  = eye(16) - K * H;
P_new = I_KH * P * I_KH' + K * R_gps * K';
P_new = (P_new + P_new') / 2;
end

%% LOCAL HELP
function dq = exp_quat(delta_theta)
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