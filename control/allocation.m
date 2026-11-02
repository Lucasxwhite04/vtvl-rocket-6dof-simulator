function [u, alloc_out] = allocation(T_cmd, alpha_cmd, beta_cmd, params)
% ALLOCATION  Combine guidance throttle and attitude gimbal commands.
%
% Sits between the guidance subsystem (provides T_cmd) and the attitude
% loop (provides alpha_cmd, beta_cmd), enforces all actuator bounds, and
% assembles the control vector u passed to the actuator models and plant.
%
% CONTROL VECTOR:
%   u = [T; alpha; beta]
%     T     — thrust magnitude (N),       clamped to [Tmin, Tmax]
%     alpha — gimbal pitch angle (rad),   clamped to ±alpha_max
%     beta  — gimbal roll  angle (rad),   clamped to ±beta_max
%
% SATURATION LOGGING:
%   alloc_out records whether any channel saturated this step, which is
%   useful for diagnosing attitude loop performance (saturation means the
%   attitude loop is asking for more than the actuator can deliver).
%
% Inputs:
%   T_cmd      scalar  commanded thrust from guidance (N)
%   alpha_cmd  scalar  commanded gimbal pitch from attitude loop (rad)
%   beta_cmd   scalar  commanded gimbal roll  from attitude loop (rad)
%   params             struct from plant_params() — uses:
%                        params.Tmin, params.Tmax
%                        params.alpha_max, params.beta_max
%
% Outputs:
%   u          (3×1)  [T; alpha; beta] — saturated actuator commands
%   alloc_out  struct  saturation flags and clipped values for logging
% L. White May 2026

% --- Clamp each channel to physical limits ---
T     = max(params.Tmin,        min(params.Tmax,       T_cmd));
alpha = max(-params.alpha_max,  min(params.alpha_max,  alpha_cmd));
beta  = max(-params.beta_max,   min(params.beta_max,   beta_cmd));

% ---- Assemble control vector ---
u = [T; alpha; beta];

% ---- Saturation flags (for logging and diagnostics) ---
alloc_out.T_sat     = (T_cmd     < params.Tmin)       || (T_cmd     > params.Tmax);
alloc_out.alpha_sat = (alpha_cmd < -params.alpha_max)  || (alpha_cmd > params.alpha_max);
alloc_out.beta_sat  = (beta_cmd  < -params.beta_max)   || (beta_cmd  > params.beta_max);
alloc_out.u         = u;

end