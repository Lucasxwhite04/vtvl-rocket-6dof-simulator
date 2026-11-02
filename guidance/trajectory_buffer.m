function [T_cmd_out, n_hat_out] = trajectory_buffer(new_trajectory, t_now, params)
% TRAJECTORY_BUFFER  The OFMS (Onboard Flight Management System) analogue.
%
% Holds the active reference trajectory. Splices in validated MPC solutions
% when they arrive. The inner attitude loop reads from this buffer at its
% own rate (100-200 Hz), independently of the MPC rate (1-2 Hz).
%
% Usage
%   trajectory_buffer('init', [], params)      Initialise from ICLOCS reference.
%                                              Call once at startup before the
%                                              first simulation step. Returns
%                                              T_cmd_out = NaN, n_hat_out = NaN.
%
%   trajectory_buffer('reset', [], params)     Alias for 'init'. Use at phase
%                                              transitions if you want a clean
%                                              buffer.
%
%   [T, n] = trajectory_buffer(sol_out, t, p) Splice sol_out into buffer, then
%                                              read the interpolated command at t.
%                                              sol_out is the struct returned by
%                                              mpc_solve.m (must have .t_grid_si,
%                                              .T_k_si, .n_hat_k fields).
%
%   [T, n] = trajectory_buffer([], t, p)      Read only at time t (no splice).
%
% Buffer implementation
%   Internally holds three parallel column vectors (buf_t, buf_T) and a
%   matrix (buf_n) kept sorted by time. All times are in seconds (SI).
%
%   On a splice, entries at or after the new trajectory's start time are
%   overwritten. This prevents stale future nodes from the previous plan
%   hanging around after a re-plan.
%
%   Reads are linear interpolations clamped to the buffer extent (no
%   extrapolation beyond the last buffer node).
%
% INPUTS:
%   new_trajectory - 'init' | 'reset' | [] | sol_out struct
%   t_now          - current mission time [s]
%   params         - vehicle params struct (must contain .reference when
%                    called as 'init'/'reset')
%
% OUTPUTS:
%   T_cmd_out      - thrust magnitude at t_now [N]
%   n_hat_out      - thrust unit direction at t_now, inertial [3x1]
%
% CALLED BY: test_step8.m, guidance_wrapper.m (Simulink)
% CALLS:     (none — self-contained with persistent state)
%
% L. White, University of Sheffield, May 2026

    % Persistent buffer: buf_t is Mx1, buf_T is Mx1, buf_n is Mx3.
    % All SI units. Sorted ascending on buf_t.
    persistent buf_t buf_T buf_n

    
    %%  INIT / RESET
    if ischar(new_trajectory) && ...
            (strcmpi(new_trajectory, 'init') || strcmpi(new_trajectory, 'reset'))

        ref = params.reference;

        % De-duplicate time axis (phase-boundary repeats break interp1)
        [t_u, ia] = unique(ref.t, 'stable');
        x_u = ref.x(ia, :);     % state rows: cols 1:7 = [r,v,m], 8:11 = q
        u_u = ref.u(ia, :);     % control rows: [T_mag, alpha_deg, beta_deg]

        M = numel(t_u);

        % Reconstruct inertial thrust unit vectors from the stored
        % control + attitude. This mirrors the logic in warm_start.m
        % ('reference' mode) so the two are consistent.
        %
        % n_body = [-sin(alpha); sin(beta); cos(alpha)*cos(beta)]
        % n_inertial = R_bi(q) * n_body
        %
        % Vectorise using the three Euler-angle–free rotation columns.
        T_mag_ref = u_u(:, 1);           % M x 1
        a_rad     = u_u(:, 2) * pi/180;  % M x 1  gimbal pitch
        b_rad     = u_u(:, 3) * pi/180;  % M x 1  gimbal roll

        n_body_mat = [-sin(a_rad), sin(b_rad), cos(a_rad).*cos(b_rad)];  % M x 3

        % Quaternion columns [q0, q1, q2, q3] from state matrix
        q0 = x_u(:, 8);  q1 = x_u(:, 9);
        q2 = x_u(:,10);  q3 = x_u(:,11);

        % Rotate each body-frame unit vector to inertial frame using R_bi.
        % R_bi rows (inertial) expressed in terms of quaternion components:
        % Row 1: [1-2(q2²+q3²), 2(q1q2-q0q3), 2(q1q3+q0q2)]
        % Row 2: [2(q1q2+q0q3), 1-2(q1²+q3²), 2(q2q3-q0q1)]
        % Row 3: [2(q1q3-q0q2), 2(q2q3+q0q1), 1-2(q1²+q2²)]
        R11 = 1 - 2*(q2.^2 + q3.^2);   R12 = 2*(q1.*q2 - q0.*q3);   R13 = 2*(q1.*q3 + q0.*q2);
        R21 = 2*(q1.*q2 + q0.*q3);      R22 = 1 - 2*(q1.^2 + q3.^2); R23 = 2*(q2.*q3 - q0.*q1);
        R31 = 2*(q1.*q3 - q0.*q2);      R32 = 2*(q2.*q3 + q0.*q1);   R33 = 1 - 2*(q1.^2 + q2.^2);

        % n_inertial = R_bi * n_body (per row, vectorised dot products)
        bx = n_body_mat(:,1);  by = n_body_mat(:,2);  bz = n_body_mat(:,3);
        n_inertial = [R11.*bx + R12.*by + R13.*bz, ...
                      R21.*bx + R22.*by + R23.*bz, ...
                      R31.*bx + R32.*by + R33.*bz];   % M x 3

        % Normalise each row (should already be unit but floating-point drift)
        row_norms = max(vecnorm(n_inertial, 2, 2), 1e-9);
        n_inertial = n_inertial ./ row_norms;

        buf_t = t_u;               % M x 1
        buf_T = T_mag_ref;         % M x 1
        buf_n = n_inertial;        % M x 3

        T_cmd_out = NaN;
        n_hat_out = NaN(3, 1);
        return;
    end

    
    %  Guard: buffer must be initialised
    if isempty(buf_t)
        error('trajectory_buffer:notInit', ...
              ['Buffer not initialised. Call trajectory_buffer(''init'', [], params) ' ...
               'before the simulation loop.']);
    end

    %  SPLICE (if a new validated MPC solution was supplied)
    if ~isempty(new_trajectory) && isstruct(new_trajectory)
        sol = new_trajectory;

        % Validate required fields
        if ~isfield(sol, 't_grid_si') || ~isfield(sol, 'T_k_si') || ~isfield(sol, 'n_hat_k')
            error('trajectory_buffer:badSol', ...
                  'sol_out must contain t_grid_si, T_k_si, and n_hat_k fields.');
        end

        t_splice = sol.t_grid_si(1);

        % Drop all buffer entries at or after the splice time.
        % Strict inequality (<) keeps the last pre-splice node as an
        % interpolation anchor, preventing a step at the join.
        keep_idx = buf_t < t_splice;
        buf_t = [buf_t(keep_idx);    sol.t_grid_si(:)];
        buf_T = [buf_T(keep_idx);    sol.T_k_si(:)];
        buf_n = [buf_n(keep_idx, :); sol.n_hat_k.'];   % n_hat_k is 3xN → (Nx3)

        % Keep sorted (sol.t_grid_si should already be monotone, but belt-and-braces)
        [buf_t, si] = sort(buf_t);
        buf_T = buf_T(si);
        buf_n = buf_n(si, :);
    end

    %  READ at t_now
    if isnan(t_now)
        T_cmd_out = NaN;
        n_hat_out = NaN(3, 1);
        return;
    end

    % Clamp to buffer extent — no extrapolation.
    % Below the buffer start: return the first entry (plan hasn't started).
    % Beyond the buffer end:  return the last entry (hold last command).
    t_q = max(buf_t(1), min(buf_t(end), t_now));

    T_cmd_out = interp1(buf_t, buf_T, t_q, 'linear');

    n_raw = interp1(buf_t, buf_n, t_q, 'linear').';   % 3x1
    if norm(n_raw) > 1e-6
        n_hat_out = n_raw / norm(n_raw);
    else
        n_hat_out = [0; 0; 1];   % degenerate fallback: point up
    end

    % Clamp thrust to physical bounds
    T_cmd_out = max(params.Tmin, min(params.Tmax, T_cmd_out));
end