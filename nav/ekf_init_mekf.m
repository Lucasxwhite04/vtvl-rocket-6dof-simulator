function [x_nom, q_ref, P] = ekf_init_mekf(x_true_0, params)
% EKF_INIT_MEKF  Initialise error-state MEKF.

% Returns:
%   x_nom (13x1)  - nominal estimates [r(3); v(3); b_g(3); b_a(3); m(1)]
%   q_ref (4x1)   - reference attitude quaternion, Hamilton scalar-first
%   P (16x16)     - error-state covariance over [δr; δv; δθ; δb_g; δb_a; δm]
% L. White, University of Sheffield, May 2026

% Nominal estimates
r_hat = x_true_0(1:3) + params.gps.pos_noise_std .* randn(3,1);
v_hat = zeros(3,1);
b_g   = zeros(3,1);
b_a   = zeros(3,1);
m_hat = params.m0 + params.ekf.mass_init_std * randn();
m_hat = max(m_hat, params.ms);
m_hat = min(m_hat, params.m0 * 1.10);

x_nom = [r_hat; v_hat; b_g; b_a; m_hat];   % (13x1)

% Reference attitude: upright on the pad 
q_ref = [1; 0; 0; 0];                       % (4x1)

% Initial error-state covariance (16x16)
sigma_r   = params.gps.pos_noise_std;
sigma_v   = [0.5; 0.5; 0.5];
sigma_att = [0.5; 0.5; 0.5] * (pi/180);     % 3x1 rad
sigma_bg  = 3 * params.imu.gyro_bias_sigma;
sigma_ba  = params.imu.accel_bias_sigma;
sigma_m   = params.ekf.mass_init_std;

% Special-case for zeroed gyro/accel bias_sigma during diagnostic runs:
% replace zero with a small floor so P doesn't have zero diagonal entries
% (which would make the filter believe biases are exactly zero forever).
sigma_bg = max(sigma_bg, 1e-3 * ones(3,1));
sigma_ba = max(sigma_ba, 1e-2 * ones(3,1));

p_diag = [sigma_r.^2;       % indices 1-3   δr
    sigma_v.^2;             % indices 4-6   δv
    sigma_att.^2;           % indices 7-9   δθ
    sigma_bg.^2;            % indices 10-12 δb_g
    sigma_ba.^2;            % indices 13-15 δb_a
    sigma_m^2];             % index   16    δm

P = diag(p_diag);


% DIAGNOSTIC: freeze biases at zero 
P(10:12, 10:12) = zeros(3,3);   % δb_g block
P(13:15, 13:15) = zeros(3,3);   % δb_a block

% Also zero any cross-correlations (should already be 0 in block-diagonal init, but be safe)
P(10:15, :)     = 0;
P(:, 10:15)     = 0;
P(10:15, 10:15) = 0;

end