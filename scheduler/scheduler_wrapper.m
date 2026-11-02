function [phase_id, t_in_phase, t_to_next] = scheduler_wrapper(x_hat, t)

% L. White, University of Sheffield, May 2026

persistent sched_state params_cached reference_cached
if isempty(params_cached)
    params_cached    = evalin('base', 'params');
    reference_cached = evalin('base', 'reference');
end
params    = params_cached;
reference = reference_cached;

% Phase transition times for 4-phase mission
t_phase = reference.phase_transition_times(2:4);

% Initialise on first call
if isempty(sched_state) || t < 0.005
    sched_state = sched_state_init(t);
end

[phase_id, t_in_phase, t_to_next, sched_state] = ...
    phase_scheduler(t, x_hat(:), t_phase, params, sched_state);

% Handle NaN t_to_next for final phase
if isnan(t_to_next)
    t_to_next = 0;
end
end