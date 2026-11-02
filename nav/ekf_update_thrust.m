function [x_hat, P] = ekf_update_thrust(x_hat, P, accel_meas, T_cmd, params)
x_hat      = x_hat(:);
accel_meas = accel_meas(:);
% EKF_UPDATE_THRUST  Pseudo-measurement update for vehicle mass.
%
% Uses the relationship:   T_cmd ≈ m · ‖f_body‖
% where f_body = accel_meas - b_a_hat (bias-corrected specific force,
% body frame). Rearranged into a scalar mass pseudo-measurement:
% L. White, University of Sheffield, May 2026

%%  BIAS-CORRECTED SPECIFIC FORCE
b_a    = x_hat(14:16);
m_hat  = x_hat(17);

% Subtract the filter's current guess for accelerometer bias, then 
% calculate the total magnitude (length) of the acceleration vector.
f_body = accel_meas - b_a;
f_norm = norm(f_body);

% A safety guard that aborts the update if acceleration is too low.
if T_cmd < 5 * params.Tmin
    return
end

% Activation gate 2: only run at low velocity to minimise drag
% contamination. At peak velocity ~30 m/s, drag ≈ 3 N which is
% small but biases the mass estimate. Below 10 m/s drag is < 0.3 N,
% negligible against typical thrust of 300-3200 N.
v_norm = norm(x_hat(4:6));
if v_norm > 10.0    % m/s threshold — tunable
    return
end

%%  PSEUDO-MEASUREMENT & OBSERVATION MATRIX
% Calculate what the mass "must be" based on the current thrust and acceleration.
z_m = T_cmd / f_norm;

% Create the Observation Matrix (H) with a '1' only in the 17th column.
H       = zeros(1, 17);
H(1,17) = 1;   % direct mass observation

%%  DYNAMIC MEASUREMENT NOISE
% Calculate how much the EKF should "distrust" this mass calculation.
sigma_accel = norm(params.ekf.accel_noise_std);   % m/s²
R_thrust    = (m_hat^2 * sigma_accel / T_cmd)^2;  % kg²

%%  KALMAN UPDATE
% Calculate the error between the virtual sensor and the EKF's current guess.
innov = z_m - m_hat;   

% Calculate the Kalman Gain to determine the optimal blending weight.
S = H * P * H' + R_thrust;   % scalar
K = P * H' / S;               % (17×1)

% Shift the 17 state estimates based on the mass error.
x_hat = x_hat + K * innov;

%%  PHYSICAL BOUNDS
% Hard-clamp the estimated mass to physical realities.
x_hat(17) = max(x_hat(17), params.ms);          % floor: dry mass
x_hat(17) = min(x_hat(17), params.m0 * 1.10);  % ceiling: 10% over nominal

%%  MATRIX SAFEGUARDS
% Force the attitude components to a length of 1.0.
x_hat(7:10) = x_hat(7:10) / norm(x_hat(7:10));

% Shrink the covariance matrix (P) using the Joseph form.
I_KH = eye(17) - K * H;
P    = I_KH * P * I_KH' + K * R_thrust * K';
P    = (P + P') / 2;   % enforce symmetry

end