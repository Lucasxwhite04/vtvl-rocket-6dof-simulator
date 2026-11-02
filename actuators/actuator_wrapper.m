function [T_act, alpha_act, beta_act] = actuator_wrapper(T_cmd, alpha_cmd, beta_cmd, t)
% L. White, University of Sheffield, May 2026
persistent valve_state act_state params_cached
if isempty(params_cached)
    params_cached = evalin('base', 'params');
end
params = params_cached;

if isempty(valve_state) || t < 0.001
    valve_state = valve_state_init(params);
    valve_state.T_act  = T_cmd;
    valve_state.T_prev = T_cmd;
    act_state = act_state_init();
    act_state.alpha_cmd_rl = alpha_cmd;
    act_state.beta_cmd_rl  = beta_cmd;
    act_state.s_alpha = [alpha_cmd; 0];
    act_state.s_beta  = [beta_cmd;  0];
end
dt = 0.005;
[T_act,    valve_state] = throttle_valve(T_cmd,   valve_state, dt, params);
[alpha_act, beta_act, act_state] = tvc_servo(alpha_cmd, beta_cmd, ...
    act_state,  dt, params);
end