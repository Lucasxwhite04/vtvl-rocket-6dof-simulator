function [phase_id, t_in_phase, t_to_next, sched_out] = phase_scheduler(t, x_hat, t_phase_nominal, params, sched_in)
% PHASE_SCHEDULER  Determines current mission phase with hysteresis.
%
% PHASE DEFINITIONS:
%   1  Ascent:   (0,0,0) → (100,-100,300) m          free duration 4-18 s
%   2  Hover:    hold (100,-100,300) m                fixed 5 s  [n_phases=4 only]
%   3  Descent:  (100,-100,300) → (200,-200,5) m      free duration ≥ 4 s
%   4  Landing:  (200,-200,5) → (200,-200,0) m        fixed 5 s
%
% Inputs:
%   t                (scalar)  current time (s)
%   x_hat            (17x1)    EKF state estimate [r;v;q;b_g;b_a;m]
%   t_phase_nominal  (3x1 or 4x1)  nominal phase end times (s)
%   params           struct    vehicle/mission params — must contain:
%                                params.n_phases  (3 or 4)
%   sched_in         struct    persistent scheduler state — init via sched_state_init(t0)
%
% Outputs:
%   phase_id    integer  1, 3, or 4 when n_phases=3; 1,2,3,4 when n_phases=4
%   t_in_phase  scalar   time elapsed since entering current phase (s)
%   t_to_next   scalar   estimated time to next transition (s); NaN in final phase
%   sched_out   struct   updated persistent state (pass as sched_in next call)
% L. White, University of Sheffield, May 2026


% --- Read hover flag from params ---
% params.n_phases = 3: ascent → descent → landing   (default)
% params.n_phases = 4: ascent → hover → descent → landing
% % HOVER_DISABLED: USE_HOVER = (params.n_phases == 4);
USE_HOVER = (params.n_phases == 4);

% --- Unpack EKF state  ---
% Scheduler only uses position (1:3) and velocity (4:6); attitude, biases,
% and mass (7:17) are not needed for phase transition logic.
r        = x_hat(1:3);
v        = x_hat(4:6);
Z        = r(3);
speed    = norm(v);
speed_xy = norm(v(1:2));

% --- Mission waypoints ---
pos_hover = [100; -100; 300];   % ascent apex / hover point (m)

% --- Scheduler parameters ---
t_dwell      = 0.3;   % hysteresis dwell time (s) - both conditions must
%                       hold continuously for this long before transition
t_hover_dur  = 4.5;   % hover fixed duration (s)           [n_phases=4 only]
t_early_gate = 2.0;   % time gate opens this many s before nominal end time

% State-gate thresholds for each transition:
pos_thresh_13 = 15.0; % ‖r_xy - r_hover_xy‖ < 15 m  (ascent apex check)
alt_thresh_13 = 260;% Z > m                     (ascent apex check)
alt_lat_thresh= 150.0; % lateral radius from apex < m
vel_thresh_13 = 4.0;  % ‖v‖ < 4 m/s        (hover-entry; n_phases=4 only)
alt_thresh_34 = 10.0;  % Z < 10 m             (descent → landing)
vel_thresh_34 = 3.0;  % ‖v_xy‖ < 3 m/s     (descent → landing)
vel_total_thresh_34 = 8.0;

% --- Map t_phase_nominal to local variables ---
if USE_HOVER
    t1_end = t_phase_nominal(1);
    t2_end = t_phase_nominal(2);
    t3_end = t_phase_nominal(3);
else
    t1_end = t_phase_nominal(1);
    t3_end = t_phase_nominal(2);
end

% --- Unpack persistent state ---
phase         = sched_in.phase;
t_entry       = sched_in.t_entry;
t_cond_on     = sched_in.t_cond_on;
cond_was_true = sched_in.cond_was_true;

t_in_phase      = t - t_entry;
dwell_satisfied = false;

% --- Evaluate transition condition for current phase ---
switch phase

    case 1  % --- Ascent - next phase ---
        if USE_HOVER
            % 1 - 2  (ascent - hover)
            % Rocket must decelerate to near-zero speed before hover entry.
            time_gate  = (t >= t1_end - t_early_gate);
            state_gate = (norm(r(1:2) - pos_hover(1:2)) < alt_lat_thresh) && ...
             (Z > alt_thresh_13) && ...
             (v(3) < 1.0);    % vz < 1 m/s — near apex or descending
            % Hard time-based fallback: prevents "stuck in Phase 1" if state gate
            % fails to fire (e.g. tracking error keeps Z just below 285 m).
            % Five-second grace period beyond the nominal phase end time.
            cond = (time_gate && state_gate) || (t > t1_end + 4.0);
            t_to_next = max(0, t1_end - t);
        else
            % 1 - 3  (ascent - descent, no hover)
            % Rocket passes through apex at speed — do NOT require low velocity.
            % Lateral XY within 20 m of apex AND altitude > 285 m prevents
            % spurious early trigger during initial climb.
            time_gate  = (t >= t1_end - t_early_gate);
            state_gate = (Z > alt_thresh_13) && (v(3) < 2.0);
            cond = (time_gate && state_gate) || (t > t1_end + 3.0);
            t_to_next = max(0, t1_end - t);
        end

   
    case 2  % --- Hover - descent ---
        % 2 - 3  time-only; hover is fixed-duration (5 s)
        cond      = (t_in_phase >= t_hover_dur);
        t_to_next = max(0, t_hover_dur - t_in_phase);

    case 3  % --- Descent → landing ---
        time_gate  = (t >= t3_end - t_early_gate);
        state_gate = (Z < alt_thresh_34) && ...
                     (norm(v) < vel_total_thresh_34);
      

        % --- absolute altitude fallback---
        % If the rocket is below 4 m, force landing regardless of speed.
        % Better to land with lateral velocity than to crash through the ground.
        alt_only_fallback = (Z < 4.0);

        % ---- update time fallbacks (fire before crash) ---
        time_fallback  = (t > t3_end + 3.0) && (Z < 30.0);   
        time_fallback2 = (t > t3_end + 12.0);                

        cond      = (time_gate && state_gate) || alt_only_fallback || ...
                    time_fallback || time_fallback2;
        t_to_next = max(0, t3_end - t);

    case 4  % --- Landing (terminal phase) ---
        cond      = false;
        t_to_next = NaN;

    otherwise
        cond      = false;
        t_to_next = NaN;
end

% --- Hysteresis de-bounce --
% Both the time gate AND state gate must hold for t_dwell seconds before
% the transition is accepted.  This prevents sensor noise or transient
% trajectory excursions from triggering early.
if cond
    if ~cond_was_true
        t_cond_on     = t;
        cond_was_true = true;
    end
    dwell_satisfied = (t - t_cond_on >= t_dwell);
else
    t_cond_on     = t;
    cond_was_true = false;
    dwell_satisfied = false;
end

% --- Apply transition ---
if dwell_satisfied
    if USE_HOVER
        next_phase = phase + 1;   % 1→2→3→4
    else
        % Skip phase 2: 1-3-4
        if phase == 1
            next_phase = 3;
        elseif phase == 3
            next_phase = 4;
        else
            next_phase = phase;
        end
    end

    if next_phase ~= phase
        phase         = next_phase;
        t_entry       = t;
        t_cond_on     = t;
        cond_was_true = false;
    end
end

% --- Pack outputs ---
phase_id   = phase;
t_in_phase = t - t_entry;

sched_out.phase         = phase;
sched_out.t_entry       = t_entry;
sched_out.t_cond_on     = t_cond_on;
sched_out.cond_was_true = cond_was_true;

end