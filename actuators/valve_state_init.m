function valve_in = valve_state_init(params)
% VALVE_STATE_INIT  Initialise throttle valve persistent state.
% Call once at simulation startup before first throttle_valve call.
% Engine starts at minimum throttle, already lit at t=0 (ascent).
% L. White, University of Sheffield, May 2026

valve_in.T_act  = params.Tmin;   % actual thrust at startup  [N]
valve_in.T_prev = params.Tmin;   % rate-limiter reference    [N]
end