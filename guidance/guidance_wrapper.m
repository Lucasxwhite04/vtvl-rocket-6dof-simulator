function [T_cmd, n_hat_cmd] = guidance_wrapper(x_hat, phase_id, t)
% GUIDANCE_WRAPPER  Mode A / Mode B dispatcher with trajectory buffer.
%
% Called extrinsically from guidance_step.m (the Simulink MATLAB Function
% block). Decides between reference-tracking (Mode A) and online MPC
% (Mode B) based on phase, manages persistent state across guidance steps,
% and uses a trajectory buffer to decouple the slow MPC planner from the
% fast inner-loop guidance rate.
%
% Architecture:
%   Phase 1, 2, 3 : Mode A — guidance_track_reference (PD over ICLOCS)
%   Phase 4       : Mode B — mpc_solve at 1 Hz; buffer interpolation between
%                            calls. Falls back to Mode A after 3 consecutive
%                            Mode B validity-check rejections (latched until
%                            the next phase transition).
%
% Persistent state (survives across guidance steps; reset on phase change):
%   prev_sol         - last accepted Mode B sol_out struct
%   rejection_count  - consecutive Mode B rejections this phase
%   fallback_active  - latched true after max_rejections; clears on
%                       phase transition
%   last_mpc_t       - time of most recent mpc_solve call
%   buf_t, buf_T, buf_n  - trajectory buffer (time, thrust mag, thrust dir)
%   last_phase_id    - phase_id from the previous call
%   first_call       - flag for one-time initialisation
%
% INPUTS (from the Simulink guidance block):
%   x_hat     - 17x1 EKF state estimate (column from Simulink)
%   phase_id  - integer 1..4 from phase scheduler
%   t         - current mission time [s]
%
% OUTPUTS:
%   T_cmd     - thrust magnitude [N]      scalar
%   n_hat_cmd - thrust unit direction     1x3 row (for Simulink logging)
%
% REQUIRES in base workspace:
%   params    - vehicle params struct, MUST INCLUDE:
%               .socp_optimizer  (compiled at startup by build_socp.m)
%               .scp.*           (SCP tolerances, weights, etc.)
%               .dt_guidance     (Mode B replan period, default 1.0 s)
%               .reference       (ICLOCS reference trajectory)
%
% CALLED BY: guidance_step.m (extrinsically)
% CALLS:     guidance_track_reference, mpc_solve
%
% L. White, University of Sheffield, May 2026

    persistent prev_sol rejection_count fallback_active last_mpc_t
    persistent buf_t buf_T buf_n last_phase_id first_call
    persistent params_cached

    %% --- One-time initialisation (and simulation reset) ---
    if isempty(first_call)
        prev_sol        = [];
        rejection_count = 0;
        fallback_active = false;
        last_mpc_t      = -inf;
        buf_t = []; buf_T = []; buf_n = [];
        last_phase_id   = 0;
        first_call      = false;
    end

    % Detect a Simulink "rerun" by t resetting to 0 — clear all state
    if t < 1e-6 && last_phase_id ~= 0
        prev_sol        = [];
        rejection_count = 0;
        fallback_active = false;
        last_mpc_t      = -inf;
        buf_t = []; buf_T = []; buf_n = [];
        last_phase_id   = 0;
    end

    %% --- Grab params and force column input ---
    if isempty(params_cached)
        params_cached = evalin('base', 'params');
    end
    params = params_cached;
    x_hat  = x_hat(:);

    %% --- Reset on phase transition ---
    if phase_id ~= last_phase_id
        prev_sol        = [];
        rejection_count = 0;
        fallback_active = false;
        last_mpc_t      = -inf;
        buf_t = []; buf_T = []; buf_n = [];
        last_phase_id   = phase_id;
    end

    %% --- Pick the active mode ---
    use_mode_b = (phase_id == 4) && ~fallback_active;

    if ~use_mode_b
        %  Mode A 
        [T_cmd, n_hat_cmd, ~] = guidance_track_reference( ...
            t, x_hat, params.reference, params);

    else
        % Mode B 
        % Trigger a replan if:
        %   (a) buffer is empty (first call in phase 4), OR
        %   (b) at least dt_guidance has elapsed since the last replan
        do_replan = isempty(buf_t) || (t - last_mpc_t >= params.dt_guidance);

        if do_replan
            target_state   = [params.r_pad; 0; 0; 0];          % pad, zero v
            sigma_bounds_s = [params.scp.sigma_land, ...
                              params.scp.sigma_land];          % fixed t_f

            [~, ~, sol_out, rej_new, fb] = mpc_solve( ...
                x_hat, t, target_state, sigma_bounds_s, ...
                prev_sol, rejection_count, params);

            last_mpc_t      = t;
            rejection_count = rej_new;

            if sol_out.valid
                prev_sol = sol_out;
                % --- Splice the new trajectory into the buffer ---
                % Drop any buffer entries from the new t_grid start onwards,
                % then concatenate and sort. This is the OFMS/PAC pattern
                % from Build Plan v2 sec 12.9.
                t_splice = sol_out.t_grid_si(1);
                keep_idx = buf_t < t_splice;
                buf_t = [buf_t(keep_idx); sol_out.t_grid_si(:)];
                buf_T = [buf_T(keep_idx); sol_out.T_k_si(:)];
                buf_n = [buf_n(keep_idx, :); sol_out.n_hat_k.'];
                [buf_t, si] = sort(buf_t);
                buf_T = buf_T(si);
                buf_n = buf_n(si, :);
            end

            if fb
                fallback_active = true;
            end
        end

        % --- Read the command from the buffer at current time ---
        if ~isempty(buf_t) && ~fallback_active
            T_cmd     = interp1(buf_t, buf_T,  t, 'linear', 'extrap');
            n_hat_cmd = interp1(buf_t, buf_n,  t, 'linear', 'extrap').';
            n_mag     = norm(n_hat_cmd);
            if n_mag > 1e-6
                n_hat_cmd = n_hat_cmd / n_mag;
            else
                n_hat_cmd = [0; 0; 1];      % defensive fallback (shouldn't fire)
            end
            % Hard clamp T_cmd to actuator limits (interp can produce
            % values microscopically outside the box).
            T_cmd = max(params.Tmin, min(params.Tmax, T_cmd));
        else
            % Either the very first mpc_solve hasn't accepted yet, OR
            % fallback got latched on this very call. Fall through to A.
            [T_cmd, n_hat_cmd, ~] = guidance_track_reference( ...
                t, x_hat, params.reference, params);
        end
    end

    %% --- Force row output for Simulink ---
    n_hat_cmd = n_hat_cmd(:)';     % 1x3
end