function r = extract_mc_result(sim_out, params, seed, scenario)
% EXTRACT_MC_RESULT  Pull per-run metrics from a Simulink sim() output.
%
% Designed to be called inside parfor — returns a small struct. Robust to
% missing logged signals: if a diagnostic signal isn't currently logged by
% vtvl_sim.slx (e.g. fallback_flag, mpc_solve_time), the field is set to NaN.
%
% INPUTS:
%   sim_out  - output of sim('vtvl_sim', ...) (or empty struct on crash)
%   params   - params struct used for this run (for r_pad, m0, thresholds)
%   seed     - run seed (recorded for traceability)
%   scenario - scenario name (string)
%
% OUTPUTS:
%   r - struct of per-run metrics. Fields:
%       .seed, .scenario, .failed, .error_msg
%       .landed                 (bool)
%       .land_x_err, .land_y_err, .land_z_err  (m, signed)
%       .land_err_radial        (m, horizontal CEP component)
%       .land_v_norm            (m/s)
%       .fuel_used              (kg)
%       .max_yaw_deg            (deg)
%       .max_pos_err_vs_truth   (m, max ‖r_hat - r_true‖)
%       .max_mass_err           (kg, max |m_hat - m_true|)
%       .max_pos_err_vs_ref     (m, vs ICLOCS reference)
%       .fallback_events        (count, NaN if signal not logged)
%       .mpc_max_solve_ms       (ms,   NaN if signal not logged)
%       .mpc_mean_iter          (—,    NaN if signal not logged)
%
% Build Plan v2 §14.4
%
% CALLED BY: run_montecarlo.m
% CALLS:     nothing
%
% REVISION HISTORY:
%   2026-XX-XX  L. White  Initial version

    % --- Initialise output with safe defaults ----------------------
    r = init_result(seed, scenario);

    try
        % --- Logged signals (mirrors run_nominal.m extraction) -----
        xt = squeeze(sim_out.x_true_log);
        if size(xt,1) < size(xt,2), xt = xt'; end

        xh = squeeze(sim_out.x_hat_log);
        if size(xh,1) < size(xh,2), xh = xh'; end

        % --- Final-state metrics -----------------------------------
        r_end_true = xt(end, 1:3)';
        v_end_true = xt(end, 4:6)';
        m_end_true = xt(end, 7);

        err_vec = r_end_true - params.r_pad;
        r.land_x_err      = err_vec(1);
        r.land_y_err      = err_vec(2);
        r.land_z_err      = err_vec(3);
        r.land_err_radial = sqrt(err_vec(1)^2 + err_vec(2)^2);
        r.land_v_norm     = norm(v_end_true);

        % "Landed" = final altitude near pad AND vehicle didn't burn through
        % its dry mass (sim that runs out of propellant and crashes through).
        % Tolerance on altitude is generous (1 m) — CEP is captured by
        % land_err_radial regardless.
        r.landed = (abs(err_vec(3)) < 1.0) && (m_end_true > params.ms * 0.99);

        r.fuel_used = params.m0 - m_end_true;

        % --- Yaw deviation (ZYX-Euler ψ from quaternion) -----------
        Q = xt(:, 8:11);                                          % [q0 q1 q2 q3]
        psi_rad = atan2(2*(Q(:,1).*Q(:,4) + Q(:,2).*Q(:,3)), ...
                        1 - 2*(Q(:,3).^2 + Q(:,4).^2));
        r.max_yaw_deg = max(abs(psi_rad)) * 180/pi;

        % --- EKF residuals -----------------------------------------
        % x_hat layout: [r(3); v(3); q(4); b_g(3); b_a(3); m(1)] = 17×1
        N_min = min(size(xt,1), size(xh,1));
        pos_err = vecnorm(xh(1:N_min, 1:3) - xt(1:N_min, 1:3), 2, 2);
        r.max_pos_err_vs_truth = max(pos_err);

        mass_err = abs(xh(1:N_min, 17) - xt(1:N_min, 7));
        r.max_mass_err = max(mass_err);

        % --- Tracking vs ICLOCS reference --------------------------
        if isfield(params, 'reference') && ~isempty(params.reference)
            ref  = params.reference;
            t_xt = linspace(0, sim_out.tout(end), size(xt,1))';
            r_ref_interp = interp1(ref.t, ref.x(:, 1:3), t_xt, 'linear', 'extrap');
            track_err = vecnorm(xt(:, 1:3) - r_ref_interp, 2, 2);
            r.max_pos_err_vs_ref = max(track_err);
        end

        % --- MPC diagnostics (optional logs) -----------------------
        % To enable, add To-Workspace blocks for these signals in vtvl_sim.slx:
        %   fallback_flag_log    (bool, 1 sample per guidance cycle)
        %   mpc_solve_time_log   (double, seconds, per guidance cycle)
        %   mpc_iter_log         (int, SCP iters per guidance cycle)
        r.fallback_events  = safe_count(sim_out, 'fallback_flag_log');
        r.mpc_max_solve_ms = safe_max(  sim_out, 'mpc_solve_time_log') * 1000;
        r.mpc_mean_iter    = safe_mean( sim_out, 'mpc_iter_log');

    catch ME
        r.failed    = true;
        r.error_msg = ME.message;
        % numeric fields already pre-filled with NaN by init_result()
    end
end

% =====================================================================
% Local helpers
% =====================================================================
function r = init_result(seed, scenario)
    r.seed     = seed;
    r.scenario = scenario;
    r.failed   = false;
    r.error_msg = '';
    r.landed   = false;

    % Pre-fill all numeric fields with NaN so downstream aggregation
    % (which does e.g. [results.land_err_radial]) doesn't crash on
    % rows that errored out.
    numeric_fields = {'land_x_err','land_y_err','land_z_err','land_err_radial', ...
                      'land_v_norm','fuel_used','max_yaw_deg', ...
                      'max_pos_err_vs_truth','max_mass_err','max_pos_err_vs_ref', ...
                      'fallback_events','mpc_max_solve_ms','mpc_mean_iter'};
    for k = 1:numel(numeric_fields)
        r.(numeric_fields{k}) = NaN;
    end
end

function n = safe_count(s, fld)
    if isstruct(s) && isfield(s, fld)
        v = squeeze(s.(fld));
        n = sum(double(v(:)) > 0.5);
    else
        n = NaN;
    end
end

function m = safe_max(s, fld)
    if isstruct(s) && isfield(s, fld)
        v = squeeze(s.(fld));
        m = max(v(:));
    else
        m = NaN;
    end
end

function m = safe_mean(s, fld)
    if isstruct(s) && isfield(s, fld)
        v = squeeze(s.(fld));
        m = mean(v(isfinite(v)));
    else
        m = NaN;
    end
end