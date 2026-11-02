%% VERIFY_EKF.m
% Step 2/3 verification: EKF estimate converges to plant truth.
%
% Test procedure:
%   1. Run plant forward for 10 s (hover) with IMU noise injected.
%   2. Run EKF alongside: predict at 200 Hz (with T_cmd), GPS update at
%      10 Hz, baro update at 50 Hz.
%   3. Check that estimation errors (position, velocity, attitude, mass)
%      are bounded and that 3-sigma covariance envelopes contain the truth.
% L. White, University of Sheffield, May 2026

clear; clc;

% Dynamically add project folders to the MATLAB path.
script_dir   = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
addpath(script_dir);
p = fullfile(project_root, 'Plant layer');
if exist(p, 'dir'), addpath(p); end
addpath(fullfile(project_root, 'Sensor Models + EKF/'));

%% Setup
params = plant_params();

% --- IMU / GPS / Baro overrides ---
% Hardcode specific sensor noise values for this test.
params.imu.rate_hz          = 200;
params.imu.accel_noise_std  = [0.01; 0.01; 0.01];
params.imu.accel_bias       = [0; 0; 0];
params.imu.accel_bias_sigma = [0.05; 0.05; 0.05];
params.imu.accel_scale_err  = [0; 0; 0];
params.imu.gyro_noise_std   = [0.001; 0.001; 0.001];
params.imu.gyro_bias        = [0; 0; 0];
params.imu.gyro_bias_sigma  = [0.005; 0.005; 0.005];
params.imu.gyro_scale_err   = [0; 0; 0];

params.gps.rate_hz          = 10;
params.gps.pos_noise_std    = [2; 2; 3];
params.gps.vel_noise_std    = [0.1; 0.1; 0.15];
params.gps.dropout_prob     = 0;

params.baro.rate_hz         = 50;
params.baro.noise_std       = 0.5;
params.baro.bias            = 0;

% --- EKF process noise ---
params.ekf.accel_noise_std       = params.imu.accel_noise_std;
params.ekf.gyro_noise_std        = params.imu.gyro_noise_std;
params.ekf.gyro_bias_rw          = [1e-4; 1e-4; 1e-4];   
params.ekf.accel_bias_rw         = [1e-3; 1e-3; 1e-3];   

% Set mass process noise and initial uncertainty.
params.ekf.mass_process_noise_std = 0.01;   
params.ekf.mass_init_std          = 0.5;    
params.ekf.thrust_vel_max         = 10.0;   

%% Simulation setup
seed      = 42;
dt_imu    = 1 / params.imu.rate_hz;   
T_sim     = 10;                       
N_steps   = round(T_sim / dt_imu);

% Calculate the ratio of IMU steps to GPS/Baro steps.
gps_ratio  = round(params.imu.rate_hz / params.gps.rate_hz);   
baro_ratio = round(params.imu.rate_hz / params.baro.rate_hz);  
wind0 = zeros(3,1);

%% Initial conditions — hovering at (100, -100, 300) m
% Define the true physical state of the rocket.
x_true  = [100; -100; 300;   
           0; 0; 0;          
           params.m0;        
           1; 0; 0; 0;       
           0; 0; 0];         
T_hover = x_true(7) * params.g0;   
u_ctrl  = [T_hover; 0; 0];         

%% EKF initialisation (17-state)
% Initialize the EKF matrices.
[x_hat, P] = ekf_init(x_true, params);

%% Preallocate logs
% Create empty arrays of the exact required size before the loop starts.
t_log       = zeros(1, N_steps);
x_true_log  = zeros(14, N_steps);  
x_hat_log   = zeros(17, N_steps);   
P_diag_log  = zeros(17, N_steps);   
innov_log   = zeros(6,  N_steps);   
qnorm_log   = zeros(1,  N_steps);
mass_cmd_log = zeros(1, N_steps);   
t = 0;

fprintf('Running EKF verification (%d steps at %d Hz) — 17-state...\n', N_steps, params.imu.rate_hz);

%% Main loop
for k = 1:N_steps
    % --- Plant truth step (midpoint Euler) ---
    %Integrate physics using the Midpoint Method (RK2).
    xdot     = plant_dynamics(x_true, u_ctrl, t, params, wind0);
    xdot_mid = plant_dynamics(x_true + 0.5*xdot*dt_imu, u_ctrl, t + 0.5*dt_imu, params, wind0);
    x_true_next = x_true + xdot_mid * dt_imu;
    
    % Update hover thrust for new (depleted) mass
    T_hover = x_true_next(7) * params.g0;
    u_ctrl  = [T_hover; 0; 0];
    
    % --- IMU measurement ---
    [accel_meas, gyro_meas] = imu_model(x_true, xdot, params, t, seed);
    
    % --- EKF predict ---
    % Run the EKF's high-speed dead-reckoning engine.
    [x_hat, P] = ekf_predict(x_hat, P, accel_meas, gyro_meas, dt_imu, params, u_ctrl(1));
    
    % --- GPS update (10 Hz) ---
    % Only run this block when 'k' reaches the gps_ratio interval.
    if mod(k-1, gps_ratio) == 0
        rng(seed + k + 6000, 'twister');
        r_meas = x_true(1:3) + params.gps.pos_noise_std .* randn(3,1);
        v_meas = x_true(4:6) + params.gps.vel_noise_std .* randn(3,1);
        innov_log(:, k) = [r_meas - x_hat(1:3); v_meas - x_hat(4:6)];
        [x_hat, P] = ekf_update_gps(x_hat, P, r_meas, v_meas, params);
    end
    
    % --- Baro update (50 Hz) ---
    if mod(k-1, baro_ratio) == 0
        rng(seed + k + 7000, 'twister');
        z_meas = x_true(3) + params.baro.bias + params.baro.noise_std * randn();
        [x_hat, P] = ekf_update_baro(x_hat, P, z_meas, params);
    end
    
    % --- Thrust pseudo-measurement ---
    % Evaluate the mass pseudo-measurement if physics conditions are safe.
    v_norm = norm(x_hat(4:6));
    if u_ctrl(1) > 5 * params.Tmin && v_norm < params.ekf.thrust_vel_max
        [x_hat, P] = ekf_update_thrust(x_hat, P, accel_meas, u_ctrl(1), params);
    end
    
    % --- Log & Advance ---
    t_log(k)         = t;
    x_true_log(:, k) = x_true;
    x_hat_log(:, k)  = x_hat;
    P_diag_log(:, k) = diag(P);
    qnorm_log(k)     = norm(x_hat(7:10));
    mass_cmd_log(k)  = u_ctrl(1);   
    x_true  = x_true_next;
    t       = t + dt_imu;
end
fprintf('Done.\n\n');

%% Extract Euler angles
% Convert quaternions into Pitch/Roll for readable plotting.
theta_hat_log  = zeros(1, N_steps);
phi_hat_log    = zeros(1, N_steps);
theta_true_log = zeros(1, N_steps);
phi_true_log   = zeros(1, N_steps);

for k = 1:N_steps
    % EKF estimate (Indices 7:10)
    q = x_hat_log(7:10, k);
    qw=q(1); qx=q(2); qy=q(3); qz=q(4);
    theta_hat_log(k) = asin( max(-1, min(1, 2*(qw*qy - qz*qx))) );
    phi_hat_log(k)   = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx^2 + qy^2));
    
    % Plant truth (Indices 8:11)
    qt = x_true_log(8:11, k);
    qw=qt(1); qx=qt(2); qy=qt(3); qz=qt(4);
    theta_true_log(k) = asin( max(-1, min(1, 2*(qw*qy - qz*qx))) );
    phi_true_log(k)   = atan2(2*(qw*qx + qy*qz), 1 - 2*(qx^2 + qy^2));
end

%% Compute errors (Ultra-Safe Base-MATLAB Math)
% Calculate the absolute difference between truth and estimate.
pos_err = x_hat_log(1:3, :) - x_true_log(1:3, :);   
vel_err = x_hat_log(4:6, :) - x_true_log(4:6, :);
att_err_theta = theta_hat_log - theta_true_log;
att_err_phi   = phi_hat_log   - phi_true_log;
mass_err = x_hat_log(17, :) - x_true_log(7, :);   

% Calculate Euclidean vector norms without relying on external toolboxes.
pos_err_norm = sqrt(sum(pos_err.^2, 1));
vel_err_norm = sqrt(sum(vel_err.^2, 1));
att_err_norm = sqrt(att_err_theta.^2 + att_err_phi.^2);

%% Pass/fail checks
fprintf('--- EKF Convergence Checks (17-state) ---\n');

% Only evaluate errors AFTER 5 seconds of flight.
idx_ss = t_log > 5;

% Calculate Root Mean Square without external toolboxes.
pos_rms  = sqrt(mean(pos_err_norm(idx_ss).^2));
vel_rms  = sqrt(mean(vel_err_norm(idx_ss).^2));
att_rms  = sqrt(mean(att_err_norm(idx_ss).^2)) * 180/pi;

mass_err_pct = mean(abs(mass_err(idx_ss)) ./ x_true_log(7, idx_ss)) * 100;
qnorm_drift  = max(abs(qnorm_log - 1));

fprintf('  Position RMS error  (steady-state t>5s): %.3f m\n',    pos_rms);
fprintf('  Velocity RMS error  (steady-state t>5s): %.4f m/s\n',  vel_rms);
fprintf('  Attitude RMS error  (steady-state t>5s): %.4f deg\n',  att_rms);
fprintf('  Mass error          (steady-state t>5s): %.3f %%\n',   mass_err_pct);
fprintf('  Quaternion norm max deviation from 1:    %.2e\n\n',    qnorm_drift);

checks.pos    = pos_rms       < 5.0;    
checks.vel    = vel_rms       < 0.5;
checks.att    = att_rms       < 2.0;
checks.mass   = mass_err_pct  < 2.0;   
checks.qnorm  = qnorm_drift   < 1e-6;

fields   = fieldnames(checks);
all_pass = true;
for i = 1:numel(fields)
    f = fields{i};
    if checks.(f)
        fprintf('  PASS: %s\n', f);
    else
        fprintf('  FAIL: %s\n', f);
        all_pass = false;
    end
end
fprintf('\n');

if all_pass
    fprintf('  ALL PASS — EKF (17-state) ready for Step 3\n\n');
else
    fprintf('  Some checks failed — tune process/measurement noise\n\n');
end

%% Plots
% Visualize the data to ensure no hidden pathological behaviors exist.
figure('Name','EKF Step 2 Verification (17-state)','Position',[60 60 1400 900]);

subplot(3,4,1);
plot(t_log, pos_err'); grid on;
xlabel('t (s)'); ylabel('Error (m)');
legend('X','Y','Z'); title('Position error');

subplot(3,4,2);
plot(t_log, vel_err'); grid on;
xlabel('t (s)'); ylabel('Error (m/s)');
legend('vx','vy','vz'); title('Velocity error');

subplot(3,4,3);
plot(t_log, [att_err_theta; att_err_phi]*180/pi); grid on;
xlabel('t (s)'); ylabel('Error (deg)');
legend('\theta error','\phi error'); title('Attitude error (Euler)');

subplot(3,4,4);
plot(t_log, mass_err); grid on; hold on;
plot(t_log, 3*sqrt(P_diag_log(17,:)), 'r--');
plot(t_log, -3*sqrt(P_diag_log(17,:)), 'r--');
xlabel('t (s)'); ylabel('Error (kg)');
legend('Mass error','3\sigma','-3\sigma'); title('Mass error with 3\sigma envelope');

subplot(3,4,5);
plot(t_log, x_true_log(1:3,:)', '--'); hold on;
plot(t_log, x_hat_log(1:3,:)', '-');
xlabel('t (s)'); ylabel('Position (m)');
legend('X true','Y true','Z true','X hat','Y hat','Z hat');
title('Position: truth vs estimate'); grid on;

subplot(3,4,6);
semilogy(t_log, sqrt(P_diag_log(1,:)), t_log, sqrt(P_diag_log(2,:)), t_log, sqrt(P_diag_log(3,:))); grid on;
xlabel('t (s)'); ylabel('1-\sigma (m)');
legend('X','Y','Z'); title('Position covariance (1-\sigma)');

subplot(3,4,7);
semilogy(t_log, sqrt(P_diag_log(11:13,:))'); grid on;
xlabel('t (s)'); ylabel('1-\sigma (rad/s)');
legend('bg_x','bg_y','bg_z'); title('Gyro bias covariance');

subplot(3,4,8);
plot(t_log, x_true_log(7,:), 'k--'); hold on;
plot(t_log, x_hat_log(17,:), 'b-'); grid on;
xlabel('t (s)'); ylabel('Mass (kg)');
legend('m_{true}','m_{hat}'); title('Mass: truth vs estimate');

subplot(3,4,9);
semilogy(t_log, abs(qnorm_log - 1) + eps); grid on;
xlabel('t (s)'); ylabel('| ||q|| - 1 |');
title('Quaternion norm deviation');

subplot(3,4,10);
semilogy(t_log, pos_err_norm, t_log, vel_err_norm, t_log, att_err_norm*180/pi); grid on;
xlabel('t (s)'); ylabel('Error norm');
legend('Pos (m)','Vel (m/s)','Att (deg)'); title('Error norms');

subplot(3,4,11);
semilogy(t_log, sqrt(P_diag_log(17,:))); grid on;
xlabel('t (s)'); ylabel('1-\sigma (kg)');
title('Mass covariance (1-\sigma)');

subplot(3,4,12);
nnz_innov = sum(any(innov_log ~= 0));
t_gps     = t_log(any(innov_log ~= 0));
innov_valid = innov_log(:, any(innov_log ~= 0));
if ~isempty(innov_valid)
    plot(t_gps, innov_valid(1:3,:)'); grid on;
    xlabel('t (s)'); ylabel('Innovation (m)');
    legend('X','Y','Z');
    title(sprintf('GPS position innovations (%d updates)', nnz_innov));
end
sgtitle('Step 2 — EKF Verification (17-state, with mass)');