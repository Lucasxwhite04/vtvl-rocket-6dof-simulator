function s = yaw_monitor_init(t0, dt)
% YAW_MONITOR_INIT  Initialise yaw monitor persistent state.
%
% Call once at mission start and pass the result as mon_in on the first
% call to yaw_monitor().
%
% Inputs:
%   t0  scalar  mission start time (s), typically 0
%   dt  scalar  monitor timestep (s) — should match attitude loop rate
% L. White May 2026

s.alarm_active = false;
s.alarm_time   = NaN;
s.alarm_psi    = NaN;
s.psi_prev     = 0;
s.t_now        = t0;
s.dt           = dt;
s.alarm_reset  = false;
end