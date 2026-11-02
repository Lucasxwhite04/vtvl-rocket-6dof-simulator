function [x_hat, P] = ekf_update_gps(x_hat, P, r_meas, v_meas, params)
% Force column vectors — defensive against Simulink shape propagation
x_hat  = x_hat(:);
r_meas = r_meas(:);
v_meas = v_meas(:);

% EKF_UPDATE_GPS  EKF measurement update with GPS position + velocity.
%
% Runs at GPS rate (10 Hz) when a valid GPS fix arrives.
% L. White, University of Sheffield, May 2026

% --- Observation matrix H (6x17) ---
% Create the mapping matrix that translates the 17-state internal 
% vector into the 6-state format that the GPS outputs.
H = zeros(6, 17);
H(1:3, 1:3) = eye(3);   % position rows
H(4:6, 4:6) = eye(3);   % velocity rows

% --- Measurement noise covariance R (6x6) ---
% Build a matrix representing how much we "distrust" the GPS data.
R_gps = diag([params.gps.pos_noise_std.^2; params.gps.vel_noise_std.^2]);

% --- Innovation (residual) ---
% Stack the 3D position and 3D velocity into a single 6x1 vector, then 
% subtract the EKF's prediction from it.
z     = [r_meas; v_meas];
z_hat = H * x_hat;        % predicted measurement (6x1)
innov = z - z_hat;        % (6x1)

% --- Innovation covariance and Kalman gain ---
% Calculate the 17x6 Kalman Gain matrix.
S = H * P * H' + R_gps;   % (6x6) Total uncertainty
K = P * H' / S;            % (17x6) Optimal blending weights

% --- Schmidt decoupling: GPS observes only position and velocity.
% Do NOT let GPS touch attitude, biases, or mass via cross-correlation.
% The 17-state EKF's additive q update is too fragile for that path.
K(7:end, :) = 0;

% --- State update ---
% Apply the mathematically weighted corrections to the current estimates.
x_hat = x_hat + K * innov;

% --- Renormalise quaternion ---
% Force the 4 attitude components back to a length of exactly 1.0.
x_hat(7:10) = x_hat(7:10) / norm(x_hat(7:10));

% --- Covariance update (Joseph form for numerical stability) ---
% Update the system's internal confidence matrix (P).
I_KH = eye(17) - K * H;
P    = I_KH * P * I_KH' + K * R_gps * K';

% Force the matrix to be perfectly symmetrical.
P    = (P + P') / 2;   % enforce symmetry

end