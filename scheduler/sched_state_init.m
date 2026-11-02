function sched_in = sched_state_init(t0)
% SCHED_STATE_INIT  Initialise phase scheduler persistent state.
% Call once at t=0 before the first phase_scheduler call.
% L. White, University of Sheffield, May 2026

sched_in.phase         = 1;      % start in ascent
sched_in.t_entry       = t0;     % time of entering phase 1
sched_in.t_cond_on     = t0;     % hysteresis timer start
sched_in.cond_was_true = false;  % no condition active yet
end