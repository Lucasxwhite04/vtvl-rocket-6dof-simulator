function [accel_meas, gyro_meas] = imu_model(x_true, xdot_true, params, t, seed)
% IMU_MODEL  Simulated IMU: 3-axis accelerometer + 3-axis gyroscope.
%
% Models body-frame specific force (accelerometer) and angular rate (gyroscope)
% with turn-on bias, white noise, and scale factor error.
% L. White, University of Sheffield, May 2026


%%  PERSISTENT BIAS — generated once per seed, not every timestep
persistent cached_seed accel_bias_rt gyro_bias_rt

p_imu = params.imu;

% Check if this is the first call or a new Monte Carlo run (new seed).
if isempty(cached_seed) || cached_seed ~= seed
    % Seed the RNG for reproducibility, adding 1000
    rng(seed + 1000, 'twister');   
    
    % Generate the permanent turn-on bias for this flight.
    accel_bias_rt = p_imu.accel_bias + p_imu.accel_bias_sigma .* randn(3, 1);
    gyro_bias_rt  = p_imu.gyro_bias  + p_imu.gyro_bias_sigma  .* randn(3, 1);
    
    cached_seed   = seed;
end

%%  UNPACK TRUE STATE  (14-state quaternion plant)
% Extract the quaternion attitude and normalize it defensively.
q0 = x_true(8);
q1 = x_true(9);
q2 = x_true(10);
q3 = x_true(11);

q_norm = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
if q_norm > eps
    q0 = q0 / q_norm;  q1 = q1 / q_norm;
    q2 = q2 / q_norm;  q3 = q3 / q_norm;
end

% Extract true physical states needed by the sensors.
p_rate     = x_true(12);   
q_rate_val = x_true(13);   
r_rate     = x_true(14);   
omega_true = [p_rate; q_rate_val; r_rate];   % (3×1) rad/s

a_cm_inertial = xdot_true(4:6);    % (3×1) m/s²
omega_dot = xdot_true(12:14);      % (3×1) rad/s²

%%  ROTATION MATRIX: body → inertial  (R_bi)  and  inertial → body (R_ib)
% Convert the quaternion into a Direction Cosine Matrix (DCM).
R_bi = [ 1 - 2*(q2^2 + q3^2),    2*(q1*q2 - q0*q3),    2*(q1*q3 + q0*q2) ;
              2*(q1*q2 + q0*q3),  1 - 2*(q1^2 + q3^2),  2*(q2*q3 - q0*q1) ;
              2*(q1*q3 - q0*q2),  2*(q2*q3 + q0*q1),    1 - 2*(q1^2 + q2^2) ];

R_ib = R_bi';   % inertial - body

%%  LEVER-ARM CORRECTION 
% Calculate the exact acceleration experienced by the physical IMU chip, 
%  which might not be located exactly at the rocket's Center of Mass (CoM).

r_imu = params.r_imu;   % (3×1) offset vector in body frame

% Step 1: Rotate CoM acceleration into the body frame
a_cm_body  = R_ib * a_cm_inertial;                         

% Step 2: Add tangential acceleration (from the rocket spinning up)
lever_angular_accel = cross(omega_dot, r_imu);             

% Step 3: Add centripetal acceleration (from the rocket holding a steady spin)
lever_centripetal   = cross(omega_true, cross(omega_true, r_imu));  

% Total IMU acceleration
a_imu_body = a_cm_body + lever_angular_accel + lever_centripetal;  


%%  TRUE SPECIFIC FORCE IN BODY FRAME
% Calculate 'Specific Force' instead of 'Coordinate Acceleration'.
g_inertial = params.grav;              % [0; 0; -9.81]
g_body     = R_ib * g_inertial;        % Rotate gravity into rocket's perspective

f_body_true = a_imu_body - g_body;     % True specific force in body frame


%%  APPLY SENSOR ERRORS
% Generate high-frequency white noise.
accel_noise = p_imu.accel_noise_std .* randn(3, 1);   
gyro_noise  = p_imu.gyro_noise_std  .* randn(3, 1);   

% Construct the final, corrupted measurement passed to the navigation computer.
accel_meas = (1 + p_imu.accel_scale_err) .* f_body_true ...
             + accel_bias_rt                             ...
             + accel_noise;                              

gyro_meas  = (1 + p_imu.gyro_scale_err)  .* omega_true ...
             + gyro_bias_rt                              ...
             + gyro_noise;                               
end