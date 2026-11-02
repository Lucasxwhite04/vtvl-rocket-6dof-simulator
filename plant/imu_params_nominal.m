function p = imu_params_nominal()
% IMU_PARAMS_NOMINAL  Representative MEMS IMU parameters for the VTVL prototype.
%
% Based on a mid-grade tactical MEMS unit (e.g. Analog Devices ADIS16488).
% Pass the returned struct into plant_params as params.imu, or use directly
% in imu_model.m.
%
% All noise values are per-sample (not per root-Hz) — they are multiplied
% by randn(3,1) inside imu_model.m to produce additive Gaussian noise.
% L. White, University of Sheffield, May 2026

% Defines the frequency at which the IMU outputs new measurements.

p.rate_hz = 200;   

% --- Accelerometer ---

% The standard deviation (1-sigma) of the high-frequency white noise.
p.accel_noise_std  = [0.01; 0.01; 0.01];    

% The expected baseline mean of the turn-on bias.
p.accel_bias       = [0;    0;    0   ];     

% The standard deviation of the permanent turn-on offset.
p.accel_bias_sigma = [0.05; 0.05; 0.05];    

% A fractional multiplier error (e.g., 0.01 = 1% error).
p.accel_scale_err  = [0;    0;    0   ];     

% --- Gyroscope ---
% White noise standard deviation for the rotational sensors.
p.gyro_noise_std   = [1e-3; 1e-3; 1e-3];   

% Nominal expected average bias of the gyro.
p.gyro_bias        = [0;    0;    0   ];    

% The standard deviation of the gyro's turn-on offset.
p.gyro_bias_sigma  = [5e-3; 5e-3; 5e-3];   

% Fractional scale error for the gyro.
p.gyro_scale_err   = [0;    0;    0   ];    

end