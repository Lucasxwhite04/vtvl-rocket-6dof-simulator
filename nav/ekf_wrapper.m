function [x_hat_out, P_diag] = ekf_wrapper(accel_meas, gyro_meas, ...
    r_meas, v_meas, gps_valid, z_meas, T_cmd, t)
T_cmd = T_cmd(1);

% L. White, University of Sheffield, May 2026

persistent x_nom q_ref P baro_counter params_cached
if isempty(params_cached)
    params_cached = evalin('base', 'params');
end
params = params_cached;
dt     = 0.005;

if isempty(x_nom) || t < dt
    x_true_0 = [0;0;0; 0;0;0; 16.88; 1;0;0;0; 0;0;0];
    [x_nom, q_ref, P] = ekf_init_mekf(x_true_0, params);
    baro_counter = 0;
end

% Predict (200 Hz)
[x_nom, q_ref, P] = ekf_predict_mekf( ...
    x_nom, q_ref, P, accel_meas(:), gyro_meas(:), dt, params, T_cmd);

% GPS update — full coupling, MEKF math
if gps_valid
    [x_nom, q_ref, P] = ekf_update_gps_mekf(x_nom, q_ref, P, ...
                         r_meas(:), v_meas(:), params);
end

% Baro update — 50 Hz: fire every 4th 200 Hz tick
baro_counter = baro_counter + 1;
if mod(baro_counter, 4) == 0
    [x_nom, q_ref, P] = ekf_update_baro_mekf(x_nom, q_ref, P, z_meas, params);
end

x_nom(13) = max(x_nom(13), params.ms);
x_nom(13) = min(x_nom(13), params.m0 * 1.10);

x_hat_full = [x_nom(1:6); q_ref; x_nom(7:12); x_nom(13)];
x_hat_out  = x_hat_full(:)';
P_diag     = [diag(P)', 0];
end