function [accel_meas, gyro_meas, r_meas, v_meas, gps_valid, z_meas] = ...
         sensors_wrapper(x_true, xdot_true, t)
% SENSORS_WRAPPER  Called as extrinsic from Simulink MATLAB Function block.
% L. White, University of Sheffield, May 2026

persistent params_cached
if isempty(params_cached)
    params_cached = evalin('base', 'params');
end
params = params_cached;

[accel_meas, gyro_meas] = imu_model(x_true, xdot_true, params, t, ...
                                     params.sim_seed);
accel_meas = accel_meas(:)';
gyro_meas  = gyro_meas(:)';
[r_meas, v_meas, gps_valid] = gps_model(x_true, params, t, params.sim_seed);
r_meas = r_meas(:)';
v_meas = v_meas(:)';
z_meas = baro_model(x_true, params, t, params.sim_seed);
end