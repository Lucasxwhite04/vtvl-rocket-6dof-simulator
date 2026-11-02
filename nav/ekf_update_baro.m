function [x_hat, P] = ekf_update_baro(x_hat, P, z_meas, params)
x_hat = x_hat(:);
% EKF_UPDATE_BARO  EKF measurement update with barometric altitude.
%
% Runs at baro rate (50 Hz). Scalar measurement: altitude only.
% L. White, University of Sheffield, May 2026

% --- Observation matrix H (1x17) ---
% Create a mapping matrix that translates the 17-state vector into what the 
% sensor should theoretically be seeing.
H_baro = zeros(1, 17);
H_baro(1, 3) = 1;   % altitude = z-component of position

% --- Measurement noise (scalar) ---
% Define the variance of the barometer.
R_baro = params.baro.noise_std^2;   % (m^2)

% --- Innovation (scalar) ---
% Calculate the difference between the actual sensor reading and what 
% the EKF *predicted* the sensor would read.
% Subtract known barometer bias before computing innovation.
% baro_model.m adds params.baro.bias to every measurement.
% Since this is a known calibration constant, we correct for it here
% rather than adding a separate bias state to the EKF.
z_corrected = z_meas - params.baro.bias;
innov = z_corrected - H_baro * x_hat;

% ---- Innovation covariance and Kalman gain ---
% Calculate the Kalman Gain (K), which is a 17x1 column vector. 
S = H_baro * P * H_baro' + R_baro;   % scalar (Total uncertainty of the measurement)
K = P * H_baro' / S;                  % (17x1) (Optimal blending factor)


% --- Schmidt decoupling: baro is altitude only.---
K(7:end) = 0;

% --- State update ---
% Shift the 17 state estimates closer to reality.
x_hat = x_hat + K * innov;

% --- Renormalise quaternion ---
% Force the 4 attitude components back to a length of exactly 1.0.
x_hat(7:10) = x_hat(7:10) / norm(x_hat(7:10));

% --- Covariance update (Joseph form for numerical stability) ---
% Shrink the uncertainty matrix (P) because we just gained new information.
I_KH = eye(17) - K * H_baro;
P    = I_KH * P * I_KH' + K * R_baro * K';

% Force the matrix to be perfectly symmetrical.
P    = (P + P') / 2;   % enforce symmetry

end