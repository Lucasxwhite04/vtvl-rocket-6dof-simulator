function [x_hat, P] = ekf_init(x_true_0, params)
% EKF_INIT  Initialise EKF state estimate and covariance from launchpad alignment.
%
% During the pre-launch alignment phase the rocket is stationary on the pad.
% We initialise position from a GPS fix, velocity as zero, attitude from
% accelerometer levelling, biases as zero (estimated in flight), and mass
% from the nominal fill value with propellant-fill uncertainty.
% L. White, University of Sheffield, May 2026


% --- Initial position: GPS fix on pad, corrupted by GPS noise ---
% Simulates taking an initial GPS reading while sitting on the pad.
r_hat = x_true_0(1:3) + params.gps.pos_noise_std .* randn(3,1);

% ---Initial velocity: stationary on pad ---
% Hardcodes the initial velocity estimate to exactly zero.
v_hat = zeros(3,1);

% --- Initial attitude: upright (identity rotation) ---
% Sets the initial attitude quaternion to [1; 0; 0; 0].
q_hat = [1; 0; 0; 0];   % [qw; qx; qy; qz]

% --- Initial biases: zero (will be estimated in flight) ---
% Sets the EKF's initial guess for the sensor turn-on biases to zero.
b_g = zeros(3,1);
b_a = zeros(3,1);

% --- Initial mass: nominal fill + random fill-uncertainty offset ----
% Guesses the starting mass by taking the theoretical nominal mass and adding 
% a random error (scaled by mass_init_std). 
m_hat = params.m0 + params.ekf.mass_init_std * randn();

% Bounds the guessed mass to physical realities.
m_hat = max(m_hat, params.ms);         
m_hat = min(m_hat, params.m0 * 1.10); 

% --- Assemble 17-state vector ---
% Stacks the individual estimates into a single 17x1 column vector.
x_hat = [r_hat; v_hat; q_hat; b_g; b_a; m_hat];

% --- Initial covariance P (17x17) ---
% Defines the initial uncertainties (variances) for all 17 states.
sigma_r   = params.gps.pos_noise_std;          % (3x1) m   - GPS position noise
sigma_v   = [0.5; 0.5; 0.5];                   % (3x1) m/s - stationary, moderate uncertainty
sigma_att = [0.01; 0.01; 0.01] * (pi/180);     % (3x1) rad - ~0.01 deg initial tilt uncertainty
sigma_bg  = 3 * params.imu.gyro_bias_sigma;    % (3x1) rad/s - 3-sigma gyro bias bound
sigma_ba  = params.imu.accel_bias_sigma;       % (3x1) m/s^2 - accel bias bound
sigma_m   = params.ekf.mass_init_std;          % (1x1) kg   - propellant fill uncertainty

% Assigns equal uncertainty to all 4 components of the quaternion.
sigma_quat = [sigma_att; sigma_att(1)];        % (4x1)

% Assembles the variances (sigma squared) into the diagonal of a 17x17 matrix.
p_diag = [sigma_r.^2;      % indices 1–3
          sigma_v.^2;      % indices 4–6
          sigma_quat.^2;   % indices 7–10
          sigma_bg.^2;     % indices 11–13
          sigma_ba.^2;     % indices 14–16
          sigma_m^2];      % index  17

P = diag(p_diag);   % (17x17)

end