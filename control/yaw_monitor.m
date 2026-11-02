function [psi_hat, alarm, mon_out] = yaw_monitor(x_hat, params, mon_in)
% YAW_MONITOR  Log yaw angle and raise alarm if it exceeds threshold.
%
% Reads the MEKF quaternion (indices 7:10 of the 17-state vector), converts
% to yaw angle ψ via quat2euler_local, and checks against the alarm
% threshold.  
%
% ALARM THRESHOLD:
%   Default: params.yaw_alarm_threshold = 10° (set in plant_params.m).
%   Alarm latches once triggered and stays high until reset via
%   mon_in.alarm_reset = true.
%
% Inputs:
%   x_hat    (17×1)  MEKF state [r; v; q; b_g; b_a; m]
%                    Quaternion at indices 7:10 (scalar-first [qw;qx;qy;qz])
%   params           struct from plant_params() — uses:
%                      params.yaw_alarm_threshold (rad)
%   mon_in   struct  persistent monitor state — init via yaw_monitor_init()
%
% Outputs:
%   psi_hat  scalar  current yaw estimate (rad)
%   alarm    bool    true if |ψ| > threshold (latches until reset)
%   mon_out  struct  updated persistent state (pass as mon_in next call)
function [psi_hat, alarm, mon_out] = yaw_monitor(x_hat, params, mon_in)
% YAW_MONITOR  Log yaw angle and raise alarm if it exceeds threshold.
%
% Reads the MEKF quaternion (indices 7:10 of the 17-state vector), converts
% to yaw angle ψ via quat2euler_local, and checks against the alarm
% threshold.  
%
% ALARM THRESHOLD:
%   Default: params.yaw_alarm_threshold = 10° (set in plant_params.m).
%   Alarm latches once triggered and stays high until reset via
%   mon_in.alarm_reset = true.
%
% Inputs:
%   x_hat    (17×1)  MEKF state [r; v; q; b_g; b_a; m]
%                    Quaternion at indices 7:10 (scalar-first [qw;qx;qy;qz])
%   params           struct from plant_params() — uses:
%                      params.yaw_alarm_threshold (rad)
%   mon_in   struct  persistent monitor state — init via yaw_monitor_init()
%
% Outputs:
%   psi_hat  scalar  current yaw estimate (rad)
%   alarm    bool    true if |ψ| > threshold (latches until reset)
%   mon_out  struct  updated persistent state (pass as mon_in next call)
% L. White May 2026

% --- Extract yaw from EKF quaternion ---
q_hat = x_hat(7:10);
[~, ~, psi_hat] = quat2euler_local(q_hat);   % discard theta, phi

% --- Threshold check ---
exceeded = abs(psi_hat) > params.yaw_alarm_threshold;

% --- Latching alarm logic ---
% Once triggered, alarm stays high until explicitly reset.
if isfield(mon_in, 'alarm_reset') && mon_in.alarm_reset
    alarm_state = false;   % reset requested
else
    alarm_state = mon_in.alarm_active || exceeded;
end

alarm = alarm_state;

% --- Log alarm event ---
if exceeded && ~mon_in.alarm_active
    % First crossing — log the time and value
    mon_out.alarm_time = mon_in.t_now;
    mon_out.alarm_psi  = psi_hat * 180/pi;   % degrees for readability
    fprintf('[YAW MONITOR] ALARM at t=%.2f s: psi=%.2f deg (threshold=%.1f deg)\n', ...
        mon_in.t_now, psi_hat*180/pi, params.yaw_alarm_threshold*180/pi);
else
    mon_out.alarm_time = mon_in.alarm_time;
    mon_out.alarm_psi  = mon_in.alarm_psi;
end

% --- Update persistent state ---
mon_out.alarm_active = alarm_state;
mon_out.psi_prev     = psi_hat;
mon_out.t_now        = mon_in.t_now + mon_in.dt;   % advance internal clock
mon_out.dt           = mon_in.dt;
mon_out.alarm_reset  = false;   % clear reset flag after one step

end


% --- Extract yaw from MEKF quaternion ---
q_hat = x_hat(7:10);
[~, ~, psi_hat] = quat2euler_local(q_hat);   % discard theta, phi

% --- Threshold check ---
exceeded = abs(psi_hat) > params.yaw_alarm_threshold;

% --- Latching alarm logic ---
% Once triggered, alarm stays high until explicitly reset.
if isfield(mon_in, 'alarm_reset') && mon_in.alarm_reset
    alarm_state = false;   % reset requested
else
    alarm_state = mon_in.alarm_active || exceeded;
end

alarm = alarm_state;

% --- Log alarm event ---
if exceeded && ~mon_in.alarm_active
    % First crossing — log the time and value
    mon_out.alarm_time = mon_in.t_now;
    mon_out.alarm_psi  = psi_hat * 180/pi;   % degrees for readability
    fprintf('[YAW MONITOR] ALARM at t=%.2f s: psi=%.2f deg (threshold=%.1f deg)\n', ...
        mon_in.t_now, psi_hat*180/pi, params.yaw_alarm_threshold*180/pi);
else
    mon_out.alarm_time = mon_in.alarm_time;
    mon_out.alarm_psi  = mon_in.alarm_psi;
end

% --- Update persistent state ---
mon_out.alarm_active = alarm_state;
mon_out.psi_prev     = psi_hat;
mon_out.t_now        = mon_in.t_now + mon_in.dt;   % advance internal clock
mon_out.dt           = mon_in.dt;
mon_out.alarm_reset  = false;   % clear reset flag after one step

end


