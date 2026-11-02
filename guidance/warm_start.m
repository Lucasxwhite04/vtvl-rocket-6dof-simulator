function [x_bar, u_bar, Gamma_bar, sigma_bar] = ...
    warm_start(mode, params, t_phase_start, t_phase_end, prev_sol)
% WARM_START  Generate the initial iterate (iterate-0) for the SCP outer loop.
%
% The SCP outer loop converges much faster from a initial guess. Two sources:
%
%   'reference' - resample the locked ICLOCS reference trajectory between
%                 t_phase_start and t_phase_end onto the SCP's N nodes.
%                 Used at the FIRST MPC call into Mode B at the start of
%                 a phase.
%
%   'shift'     - shift the previous accepted MPC solution forward by
%                 params.dt_guidance and re-pad the tail to maintain N
%                 nodes. Used for subsequent MPC calls within the same
%                 phase.
%
% In both cases the output is in non-dimensional units, ready to hand to
% scp_iteration.m as (x_bar, u_bar, Gamma_bar, sigma_bar).
%
% --- Inputs--
%   mode           - 'reference' or 'shift'
%   params         - vehicle params struct (must contain .reference,
%                    .scp.N, U_L, U_V, U_M, U_F, U_T, dt_guidance)
%   t_phase_start  - reference time at which this phase begins  [s]
%                    (used only in 'reference' mode; pass [] for 'shift')
%   t_phase_end    - reference time at which this phase ends    [s]
%                    (used only in 'reference' mode; pass [] for 'shift')
%   prev_sol       - previous MPC solution struct (used only when
%                    mode = 'shift'). Must contain:
%                      .x_k      7 x N  (non-dim)
%                      .u_k      3 x N  (non-dim)
%                      .Gamma_k  1 x N  (non-dim)
%                      .sigma    scalar (non-dim)
%                    Pass [] if mode = 'reference'.
%
% -- Outputs (all non-dimensional) ---
%   x_bar      - state iterate     7 x N
%   u_bar      - thrust vector     3 x N
%   Gamma_bar  - LCvx slack        1 x N
%   sigma_bar  - time dilation     scalar
%                NOTE: mpc_solve.m overrides sigma_bar with
%                mean(sigma_bounds_s)/U_T after this call.
%                sigma_bar returned here is "best guess" for direct
%                callers (e.g. test scripts).
%
% CALLED BY: mpc_solve.m (and the SCP-loop test scripts)
% CALLS:     nondim_state.m
% L. White, University of Sheffield, May 2026

    N = params.scp.N;

    switch lower(mode)
    %%  'reference' mode - sample from ICLOCS reference trajectory
   
        case 'reference'
            ref = params.reference;
            sigma_dim = t_phase_end - t_phase_start;
            sigma_bar = sigma_dim / params.U_T;        % non-dim

            % De-duplicate phase-boundary repeats (interp1 needs unique x)
            [t_unique, ia] = unique(ref.t, 'stable');
            x_unique = ref.x(ia, :);
            u_unique = ref.u(ia, :);

            % Resample to N nodes
            t_grid   = linspace(t_phase_start, t_phase_end, N).';
            x_phase   = interp1(t_unique, x_unique(:, 1:7), t_grid, 'linear', 'extrap');
            u_phase   = interp1(t_unique, u_unique,         t_grid, 'linear', 'extrap');

            % Pre-interpolate the quaternion history for thrust-vector
            % reconstruction (reference stores u as [T, alpha, beta]).
            q_all     = x_unique(:, 8:11);
            q_at_grid = interp1(t_unique, q_all, t_grid, 'linear', 'extrap');

            % Reconstruct the inertial thrust vector at each node:
            %   body_frame n = [-sin(alpha); sin(beta); cos(alpha)cos(beta)]
            %   inertial = R_bi(q) * (T * n_body)
            T_vec_si = zeros(3, N);
            Gamma_si = zeros(1, N);
            for k = 1:N
                T_mag = u_phase(k, 1);
                a_rad = u_phase(k, 2) * pi/180;
                b_rad = u_phase(k, 3) * pi/180;
                n_body = [-sin(a_rad); sin(b_rad); cos(a_rad)*cos(b_rad)];

                q0 = q_at_grid(k,1); q1 = q_at_grid(k,2);
                q2 = q_at_grid(k,3); q3 = q_at_grid(k,4);
                R_bi = [1-2*(q2^2+q3^2),  2*(q1*q2-q0*q3),  2*(q1*q3+q0*q2);
                        2*(q1*q2+q0*q3),  1-2*(q1^2+q3^2),  2*(q2*q3-q0*q1);
                        2*(q1*q3-q0*q2),  2*(q2*q3+q0*q1),  1-2*(q1^2+q2^2)];

                T_vec_si(:,k) = R_bi * (T_mag * n_body);
                Gamma_si(k)   = T_mag;       % LCvx tight on the optimum
            end

            % Convert to non-dim
            x_bar     = nondim_state(x_phase.', params, 'to_nd');     % 7 x N
            u_bar     = T_vec_si / params.U_F;                         % 3 x N
            Gamma_bar = Gamma_si / params.U_F;                         % 1 x N

    %%    'shift' mode - shift previous solution forward by dt_guidance
        case 'shift'
            % Concept: dt_guidance seconds have elapsed since the previous
            % MPC call. Some leading nodes of prev_sol are now in the past
            % (already executed). Drop them and pad the tail with a
            % copy of the last node so the array stays N columns wide.
            %
            % The shifted iterate is still non-dimensional and is handed
            % directly to scp_iteration.m. The SCP will reshape it toward
            % the new initial condition x0_nd and new sigma_bounds.

            if isempty(prev_sol)
                error('warm_start:noPrev', ...
                      'mode=''shift'' requires a non-empty prev_sol struct.');
            end

            % --- Unpack previous solution (all non-dim) ---
            x_prev     = prev_sol.x_k;      % 7 x N
            u_prev     = prev_sol.u_k;      % 3 x N
            Gamma_prev = prev_sol.Gamma_k;  % 1 x N
            sigma_prev = prev_sol.sigma;    % scalar (non-dim)

            sigma_si_prev = sigma_prev * params.U_T;   % previous horizon [s]

            % --- Compute how many nodes to drop ---
            % The previous trajectory covered N nodes over sigma_si_prev
            % seconds, so adjacent nodes are spaced:
            %   node_dt = sigma_si_prev / (N - 1)   [s]
            %
            % In dt_guidance seconds we've consumed approximately:
            %   n_shift = floor(dt_guidance / node_dt)
            %           = floor(dt_guidance * (N-1) / sigma_si_prev)
            %
            % Guard rails:
            %   min 1  - always make at least one step of progress
            %   max N-2 - always keep at least 2 original nodes so
            %             the SCP has a meaningful initial trajectory
            node_dt_s = sigma_si_prev / (N - 1);
            n_shift   = floor(params.dt_guidance / node_dt_s);
            n_shift   = max(1, min(n_shift, N - 2));

            % --- Shift and pad ---
            % Keep nodes [n_shift+1 .. N], repeat the last node n_shift times.
            kept = N - n_shift;          % number of kept columns (>= 2)

            x_bar     = [x_prev(:, n_shift+1:end),    repmat(x_prev(:,end),    1, n_shift)];
            u_bar     = [u_prev(:, n_shift+1:end),    repmat(u_prev(:,end),    1, n_shift)];
            Gamma_bar = [Gamma_prev(n_shift+1:end),   repmat(Gamma_prev(end),  1, n_shift)];

            % Sanity: output must be exactly 7xN, 3xN, 1xN
            assert(size(x_bar,2) == N && size(u_bar,2) == N && numel(Gamma_bar) == N, ...
                'warm_start:shiftSizeMismatch', ...
                'Shift produced wrong column count (kept=%d, n_shift=%d, N=%d).', ...
                kept, n_shift, N);

            % --- sigma_bar: best-guess remaining horizon ---
            % The caller (mpc_solve.m) overrides this with
            %   mean(sigma_bounds_s) / U_T
            % but we return a sensible value for direct test-script callers.
            sigma_remaining_s = max(sigma_si_prev - params.dt_guidance, ...
                                    0.1 * params.U_T);
            sigma_bar = sigma_remaining_s / params.U_T;

        otherwise
            error('warm_start:badMode', ...
                  'mode must be ''reference'' or ''shift'' (got "%s").', mode);
    end
end