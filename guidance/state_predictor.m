function x_predicted_si = state_predictor(x7_si, t_now, t_apply, prev_sol, params)
% STATE_PREDICTOR  G-FOLD-pattern state predictor for planning latency.
%
% The MPC takes ~30-100 ms to solve. By the time we hand a trajectory
% back to the attitude loop, the vehicle has moved. If we plan from the
% state at t_now, the trajectory will have a step discontinuity at the
% splice point. The predictor eliminates this by projecting the EKF
% estimate forward to t_apply = t_now + t_plan_delay using the most
% recent thrust history.
%
% Two regimes:
%   - First MPC call (prev_sol = []):  gravity-only propagation.
%     This is conservative - if the vehicle is in fact at trim hover,
%     we underestimate altitude slightly; nothing catastrophic.
%   - Subsequent calls (prev_sol present):  integrate the 3-DoF
%     translational dynamics using the thrust commanded by the previous
%     solution at t_now (interpolated linearly within the prev_sol grid).
%     Drag is included via transl_dynamics.
%
% --- Inputs ---
%   x7_si      - 7x1 state in SI: [r(3); v(3); m]
%                  positions [m], velocities [m/s], mass [kg]
%   t_now      - time at which the EKF estimate was formed [s]
%   t_apply    - time at which the new trajectory will start [s]
%   prev_sol   - previous MPC solution struct (or [] for first call)
%                  needs:  .u_k        (3 x N, non-dim thrust vector)
%                          .Gamma_k    (1 x N, non-dim LCvx slack)
%                          .t_grid_si  (1 x N, time stamps in seconds)
%   params     - vehicle params (uses U_F, grav, dispersions, etc.)
%
% --- Outputs ---
%   x_predicted_si - 7x1 SI state at t_apply
%
% CALLED BY: mpc_solve.m
% CALLS:     transl_dynamics.m
% L. White, University of Sheffield, May 2026

dt = t_apply - t_now;
if dt < 0
    error('state_predictor:badTime', ...
        't_apply (%.3f) must be >= t_now (%.3f).', t_apply, t_now);
end
if dt < 1e-9
    % Zero-delay case: nothing to predict
    x_predicted_si = x7_si;
    return;
end

if isempty(prev_sol) || ~isfield(prev_sol, 'u_k')
    % --- First MPC call: gravity-only propagation ---
    % r(t+dt) = r + v*dt + 0.5*g*dt^2
    % v(t+dt) = v + g*dt
    % m(t+dt) = m  (no thrust, no mass change)
    r0 = x7_si(1:3);  v0 = x7_si(4:6);  m0 = x7_si(7);
    r_p = r0 + v0*dt + 0.5*params.grav*dt^2;
    v_p = v0 + params.grav*dt;
    m_p = m0;
    x_predicted_si = [r_p; v_p; m_p];
    return;
end

% --- Subsequent calls: integrate with previously commanded thrust ---
% Interpolate prev_sol's thrust vector and LCvx slack at t_now.
% prev_sol is in non-dim units; convert to SI for transl_dynamics.
t_grid = prev_sol.t_grid_si(:);                         % column for interp1
u_nd_at_t = interp1(t_grid, prev_sol.u_k.',     t_now, 'linear', 'extrap').';
G_nd_at_t = interp1(t_grid, prev_sol.Gamma_k.', t_now, 'linear', 'extrap');

T_vec_si = u_nd_at_t * params.U_F;     % to Newtons
Gamma_si = G_nd_at_t * params.U_F;

% Integrate the 3-DoF model from t_now to t_apply.
% Loose tolerances are fine - this is a 100 ms propagation.
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
odefun = @(~, x) transl_dynamics(x, T_vec_si, Gamma_si, params);
[~, Xp] = ode45(odefun, [t_now, t_apply], x7_si, ode_opts);

x_predicted_si = Xp(end, :).';
end