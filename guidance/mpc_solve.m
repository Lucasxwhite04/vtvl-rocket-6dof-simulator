function [T_cmd, n_hat_cmd, sol_out, rejection_count_new, fallback_flag] = ...
    mpc_solve(x_hat, t_now, target_state, sigma_bounds_s, ...
              prev_sol, rejection_count, params)
% MPC_SOLVE  Top-level Mode B entry. One MPC call -> one command + trajectory.
%
% Called by the Guidance Subsystem at the guidance rate (1-2 Hz). Wires
% together the entire Mode B pipeline:
%
%   1. Slice the 7-state translational [r; v; m] from the 17-state MEKF estimate
%   2. State predictor: project forward by params.t_plan_delay
%   3. Warm start: from reference (first call in phase) or shift (within phase)
%   4. SCP outer loop until cost stabilises or i_max iterations hit
%   5. Re-dimensionalise the converged trajectory
%   6. Assemble sol_out struct (full schema, see below)
%   7. Extract first-step command via extract_command.m
%   8. Validity gate via validity_check.m -> sets fallback_flag if needed
%
% --- Inputs ---
%   x_hat            17 x 1   MEKF state estimate [SI]
%                              (1:3)=r [m], (4:6)=v [m/s], (7:10)=q,
%                              (11:13)=b_g, (14:16)=b_a, (17)=mass [kg]
%   t_now            scalar   current mission time [s]
%   target_state     6 x 1    [r_target(3); v_target(3)] in SI
%                              (mass left free for the SOCP to minimise)
%   sigma_bounds_s   1 x 2    [sigma_low, sigma_high] in seconds
%                              landing: [3, 3]   descent: [5, 30]
%   prev_sol         struct   previous accepted solution OR []
%   rejection_count  integer  consecutive Mode B rejections (this phase)
%   params           struct   vehicle params (must contain .socp_optimizer
%                              -- compiled at startup by build_socp.m)
%
% --- Outputs ---
%   T_cmd                scalar   first-step thrust magnitude [N]
%   n_hat_cmd            3 x 1    first-step thrust direction (inertial)
%   sol_out              struct   full solution (see field list below)
%   rejection_count_new  integer  updated rejection counter
%   fallback_flag        bool     true if max_rejections reached
%                                 (caller should switch to Mode A)
%
% --- sol_out schema ---
%   .x_k            7 x N   (non-dim)  state trajectory
%   .u_k            3 x N   (non-dim)  thrust vector
%   .Gamma_k        1 x N   (non-dim)  LCvx slack
%   .sigma          scalar  (non-dim)  time dilation
%   .x_target       6 x 1   (non-dim)  target supplied to SCP
%   .t_grid_si      1 x N   [s]        dimensional time stamps
%   .T_k_si         1 x N   [N]        thrust magnitudes
%   .n_hat_k        3 x N   (unit)     thrust directions
%   .lcvx_slack_N   scalar  [N]        max LCvx slack across nodes
%   .nu_norm        scalar             max ||nu_k|| from final iter
%   .exit_code      scalar             YALMIP code from final iter
%   .iter           scalar             # SCP iterations actually used
%   .iter_history   struct array       per-iter diagnostics
%   .valid          bool               from validity_check
%   .fault_code     scalar             from validity_check
%   .t_solve_total  scalar  [s]        wall-clock for the whole call
%
% CALLED BY: Simulink guidance dispatcher (or test scripts)
% CALLS:     state_predictor, warm_start, scp_iteration,
%            extract_command, validity_check, nondim_state
% L. White, University of Sheffield, May 2026

    t_start = tic;

    %% Step 1 - slice the 7-state translational from the 17-state EKF
    x7_si = [x_hat(1:6); x_hat(17)];           % [r; v; m] in SI

    %% Step 2 - state predictor: project forward by t_plan_delay
    t_apply = t_now + params.t_plan_delay;
    x0_si = state_predictor(x7_si, t_now, t_apply, prev_sol, params);
    x0_nd = nondim_state(x0_si, params, 'to_nd');           % 7 x 1, non-dim

    %% Step 3 - determine the phase from t_apply (selects warm-start segment)
    phase_id = identify_phase(t_apply, params.reference.phase_transition_times);

    %% Step 4 - warm start
    %   - first call in this phase, OR a different phase than prev_sol, OR
    %     no prev_sol available  ->  warm-start from reference
    %   - else                       ->  shift (deferred to Step 8)
    need_reference_warmstart = isempty(prev_sol) ...
        || ~isfield(prev_sol, 'phase_id') ...
        || prev_sol.phase_id ~= phase_id;

   
    if need_reference_warmstart
        phase_end_t = params.reference.phase_transition_times(phase_id + 1);
        ws_start_t  = t_apply;
        ws_end_t    = min(phase_end_t, t_apply + sigma_bounds_s(2));

        if ws_end_t <= ws_start_t + 1e-3
            ws_end_t = ws_start_t + max(sigma_bounds_s(2), 0.5);
        end

        [x_bar, u_bar, Gamma_bar, sigma_bar] = warm_start( ...
            'reference', params, ws_start_t, ws_end_t, []);
    else
        % 'shift' mode (within-phase re-plans)
        [x_bar, u_bar, Gamma_bar, sigma_bar] = warm_start( ...
            'shift', params, [], [], prev_sol);
    end

    % Sigma warm-start: applies to BOTH reference and shift modes.
    % Use the previous converged sigma if available; otherwise cold-start
    % at the midpoint of the bounds. Clamp to the bounds either way.
    sigma_bounds_nd = sigma_bounds_s / params.U_T;
    if ~isempty(prev_sol) && isfield(prev_sol, 'sigma') && ~isempty(prev_sol.sigma)
        sigma_bar = prev_sol.sigma;
    else
        sigma_bar = mean(sigma_bounds_s) / params.U_T;
    end
    sigma_bar = max(sigma_bounds_nd(1), min(sigma_bounds_nd(2), sigma_bar));

    %% Step 5 - boundary conditions for the SOCP
    % x0_nd is the predicted (post-100ms) initial condition.
    % x_target is target_state non-dimensionalised.
    r_target_nd = target_state(1:3) / params.U_L;
    v_target_nd = target_state(4:6) / params.U_V;
    x_target_nd = [r_target_nd; v_target_nd];               % 6 x 1

    sigma_low_nd  = sigma_bounds_s(1) / params.U_T;
    sigma_high_nd = sigma_bounds_s(2) / params.U_T;

    % Time-penalty weight: 0 for landing (fixed t_f), 1 for descent
    if abs(sigma_bounds_s(1) - sigma_bounds_s(2)) < 1e-6
        w_t_val = params.scp.w_t_land;                       % fixed sigma
    else
        w_t_val = params.scp.w_t_desc;                       % free sigma
    end

    %% Step 6 - SCP outer loop (cost-based termination)
    i_max   = params.scp.i_max;
    nu_tol  = params.scp.nu_tol;
    if isfield(params.scp, 'cost_tol')
        cost_tol = params.scp.cost_tol;
    else
        cost_tol = 1e-4;                                     % 0.1 g fuel
    end

    iter_hist(i_max) = struct('nu_norm',[],'delta_norm',[], ...
                              'lcvx_slack',[],'fuel_kg',[], ...
                              't_lin',[],'t_solve',[],'exit_code',[]);
    prev_fuel_kg = inf;
    converged = false;
    iters_used = 0;

    for i = 1:i_max
        [x_new, u_new, Gamma_new, sigma_new, info] = scp_iteration( ...
            x_bar, u_bar, Gamma_bar, sigma_bar, ...
            x0_nd, x_target_nd, sigma_low_nd, sigma_high_nd, w_t_val, ...
            params.socp_optimizer, params);

        fuel_kg = (1 - x_new(7,end)) * params.U_M;           % m0 - m_f

        iter_hist(i).nu_norm    = info.nu_norm;
        iter_hist(i).delta_norm = info.delta_norm;
        iter_hist(i).lcvx_slack = info.lcvx_slack;
        iter_hist(i).fuel_kg    = fuel_kg;
        iter_hist(i).t_lin      = info.t_lin;
        iter_hist(i).t_solve    = info.t_solve;
        iter_hist(i).exit_code  = info.exit_code;

        iters_used = i;

        % Bail if the solver couldn't return a valid solution at all
        if info.exit_code ~= 0
            break;
        end

        % Update the iterate
        x_bar     = x_new;
        u_bar     = u_new;
        Gamma_bar = Gamma_new;
        sigma_bar = sigma_new;

        % Cost-based termination (primary) gated on nu being small
        if i > 1 && info.nu_norm < nu_tol
            cost_change = abs(fuel_kg - prev_fuel_kg);
            if cost_change < cost_tol
                converged = true;
                break;
            end
        end
        prev_fuel_kg = fuel_kg;
    end

    iter_hist = iter_hist(1:iters_used);                     % trim unused

    %% Step 7 - re-dimensionalise and assemble sol_out
    % Build dimensional time grid: nodes equally spaced between t_apply
    % and t_apply + sigma_si.
    sigma_si = sigma_new * params.U_T;
    t_grid_si = linspace(t_apply, t_apply + sigma_si, params.scp.N);

    % Thrust magnitudes and unit directions per node (SI)
    u_si = u_new * params.U_F;                                % 3 x N, [N]
    T_k_si = vecnorm(u_si);                                   % 1 x N, [N]
    n_hat_k = u_si ./ max(T_k_si, 1e-9);                      % 3 x N, unit

    sol_out = struct( ...
        'phase_id',      phase_id, ...
        'x_k',           x_new, ...
        'u_k',           u_new, ...
        'Gamma_k',       Gamma_new, ...
        'sigma',         sigma_new, ...
        'x_target',      x_target_nd, ...
        't_grid_si',     t_grid_si, ...
        'T_k_si',        T_k_si, ...
        'n_hat_k',       n_hat_k, ...
        'lcvx_slack_N',  info.lcvx_slack, ...
        'nu_norm',       info.nu_norm, ...
        'exit_code',     info.exit_code, ...
        'iter',          iters_used, ...
        'iter_history',  iter_hist, ...
        'converged',     converged, ...
        'valid',         false, ...                           % set by check
        'fault_code',    -1, ...                              % set by check
        't_solve_total', 0);                                  % set below

    %% Step 8 - extract the first-step command
    [T_cmd, n_hat_cmd, ~] = extract_command( ...
        u_new(:,1), Gamma_new(1), params);

    %% Step 9 - validity gate
    [is_valid, fault_code, fallback_flag, rejection_count_new] = ...
        validity_check(sol_out, rejection_count, params);

    sol_out.valid         = is_valid;
    sol_out.fault_code    = fault_code;
    sol_out.t_solve_total = toc(t_start);
end


%%  LOCAL - phase identification from current time
function phase_id = identify_phase(t, transitions)
% Returns the integer phase index containing time t.
% Uses the half-open interval [t_k, t_{k+1}).  Clamps at edges.
    if t <= transitions(1)
        phase_id = 1;
        return;
    end
    if t >= transitions(end)
        phase_id = numel(transitions) - 1;
        return;
    end
    phase_id = find(t >= transitions(1:end-1) & t < transitions(2:end), 1);
end