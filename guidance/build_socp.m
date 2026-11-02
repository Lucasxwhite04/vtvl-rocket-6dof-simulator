function optimizer_obj = build_socp(params)
% BUILD_SOCP  Compile the parametric SOCP for Mode B once at startup.
%
% Returns a YALMIP optimizer object that mpc_solve.m calls with numeric
% values per SCP iteration. Compiling once and re-using the handle is
% what keeps the per-call solve under 50 ms.
%
% --- Decision variables (non-dim) ---
%   x       7 x N      state trajectory [r; v; m]
%   u       3 x N      thrust vector T_vec
%   Gamma   1 x N      LCvx thrust-magnitude slack  (||T_vec|| <= Gamma)
%   sigma   1 x 1      free-final-time dilation (sigma = t_f / U_T)
%   nu      7 x (N-1)  virtual control (dynamics defect, soft)
%   s_delta 1 x N      trust-region slack per node (soft)
%
% --- Runtime parameters (in this order) ---
%   Abar         7 x 7 x (N-1)   discrete state Jacobian
%   Bbar         7 x 4 x (N-1)   discrete control Jacobian (start)
%   Cbar         7 x 4 x (N-1)   discrete control Jacobian (end, FOH)
%   Sigbar       7 x (N-1)       sigma column
%   zbar         7 x (N-1)       affine remainder
%   x0           7 x 1           initial condition (predicted state)
%   x_target     6 x 1           [r_target; v_target] (mass is FREE)
%   x_prev       7 x N           trust-region centre (state iterate)
%   u_prev       3 x N           trust-region centre (control iterate)
%   Gamma_prev   1 x N           trust-region centre (slack iterate)
%   sigma_low    1 x 1           lower bound on sigma
%   sigma_high   1 x 1           upper bound on sigma
%   w_t          1 x 1           time-penalty weight (0 landing, 1 descent)
%
% --- Outputs (returned by optimizer_obj(...) call) ---
%   {x, u, Gamma, sigma, nu, s_delta}
%
% --- Constants baked in at compile time ----
%   T_min_nd, T_max_nd      thrust bounds (non-dim)
%   m_dry_nd                dry-mass floor (non-dim)
%   r_pad_nd                landing-pad position (non-dim)
%   gamma_gs                glideslope half-angle (rad)
%   delta_max               trust-region radius (non-dim)
%   w_nu, w_delta           virtual-control & trust-region weights
%
% --- Constraints ---
%   (1) Initial:        x(:,1) == x0
%   (2) Terminal:       x(1:6, N) == x_target           (mass free)
%   (3) Dynamics FOH:   x(:,k+1) == Abar*x_k + Bbar*u_aug_k + Cbar*u_aug_{k+1} + Sigbar*sigma + zbar + nu_k
%                       where u_aug = [T_vec; Gamma]
%   (4) LCvx SOC:       ||u(:,k)|| <= Gamma(k)
%   (5) LCvx box:       T_min_nd <= Gamma(k) <= T_max_nd
%   (6) Glideslope SOC: ||x(1:2,k) - r_pad(1:2)|| <= tan(gamma_gs) * (x(3,k) - r_pad(3))
%   (7) Subsurface:     x(3,k) >= r_pad(3)
%   (8) Mass floor:     x(7,k) >= m_dry_nd
%   (9) Sigma bounds:   sigma_low <= sigma <= sigma_high
%  (10) Trust region:   ||[x_k - x_prev_k; u_k - u_prev_k; Gamma_k - Gamma_prev_k]||
%                       <= delta_max + s_delta(k);   s_delta(k) >= 0
% --- Cost ---
%   J = -x(7,N)                  (maximise final mass = minimise propellant)
%       + w_nu * ||nu||_1
%       + w_delta * sum(s_delta)
%       + w_t * sigma
% Why minimise -m_f (and not sum(Gamma)) for fuel:
%   The true fuel cost is sigma * sum(Gamma) * dtau / c, which is BILINEAR
%   in sigma and Gamma (not convex). Minimising -m_f is convex (m_f is a
%   state variable) and exactly equivalent because the dynamics conserve
%   mass: m_f = m_0 - propellant_burned.
%
% CALLED BY: mpc_solve.m (once at startup, stored in params.socp_optimizer)
% CALLS:     none (uses YALMIP + ECOS)
% L. White, University of Sheffield, May 2026



    %  Step 1 - Extract sizes and bake-in constants from params
    
    N = params.scp.N;

    T_min_nd = params.Tmin / params.U_F;        % LCvx box lower
    T_max_nd = params.Tmax / params.U_F;        % LCvx box upper
    m_dry_nd = params.ms   / params.U_M;        % dry-mass floor
    r_pad_nd = params.r_pad / params.U_L;       % landing pad (non-dim)
    gamma_gs = params.gamma_gs;                 % radians
    tan_gs   = tan(gamma_gs);                   % bake the tan once

    delta_max = params.scp.delta_max;           % trust-region radius
    w_nu      = params.scp.w_nu;                % virtual-control penalty
    w_delta   = params.scp.w_delta;             % trust-region slack penalty



    %  Step 2 - Clear YALMIP and declare decision variables
    yalmip('clear');

    x       = sdpvar(7, N, 'full');             % state trajectory
    u       = sdpvar(3, N, 'full');             % thrust vector T_vec
    Gamma   = sdpvar(1, N, 'full');             % LCvx slack
    sigma   = sdpvar(1, 1);                     % time dilation
    nu      = sdpvar(7, N-1, 'full');           % virtual control
    s_delta = sdpvar(1, N, 'full');             % trust-region slack


    %  Step 3 - Declare runtime-parameter sdpvars
    Abar       = sdpvar(7, 7, N-1, 'full');
    Bbar       = sdpvar(7, 4, N-1, 'full');
    Cbar       = sdpvar(7, 4, N-1, 'full');
    Sigbar     = sdpvar(7, N-1, 'full');
    zbar       = sdpvar(7, N-1, 'full');

    x0         = sdpvar(7, 1);
    x_target   = sdpvar(6, 1);                  % [r_target; v_target]
    x_prev     = sdpvar(7, N, 'full');
    u_prev     = sdpvar(3, N, 'full');
    Gamma_prev = sdpvar(1, N, 'full');
    sigma_low  = sdpvar(1, 1);
    sigma_high = sdpvar(1, 1);
    w_t        = sdpvar(1, 1);


    %  Step 4 - Build constraints 
    C = [];

    % --- (1) Initial condition ---
    C = [C; x(:,1) == x0];

    % --- (2) Terminal: pin r and v, leave m free ---
    C = [C; x(1:6, N) == x_target];

    % --- (3) Dynamics (FOH on augmented control u_aug = [T_vec; Gamma]) --
    for k = 1:N-1
        u_aug_k   = [u(:,k);   Gamma(k)];
        u_aug_kp1 = [u(:,k+1); Gamma(k+1)];
        C = [C; x(:,k+1) == Abar(:,:,k) * x(:,k) ...
                          + Bbar(:,:,k) * u_aug_k ...
                          + Cbar(:,:,k) * u_aug_kp1 ...
                          + Sigbar(:,k) * sigma ...
                          + zbar(:,k) ...
                          + nu(:,k)];
    end

    % --- (4) LCvx SOC: ||T_vec_k|| <= Gamma_k ---
    %
    % --- (5) LCvx box: T_min_nd <= Gamma_k <= T_max_nd ---
    for k = 1:N
        C = [C; norm(u(:,k)) <= Gamma(k)];
        C = [C; T_min_nd <= Gamma(k) <= T_max_nd];
    end

    % --- (6) Glideslope SOC ---
    %   ||r_xy - r_pad_xy|| <= tan(gamma_gs) * (r_z - r_pad_z)
    for k = 1:N
        C = [C; norm(x(1:2,k) - r_pad_nd(1:2)) ...
                  <= tan_gs * (x(3,k) - r_pad_nd(3))];
    end

    % --- (7) Subsurface ---
    for k = 1:N
        C = [C; x(3,k) >= r_pad_nd(3)];
    end

    % --- (8) Mass floor ---
    % Without this the solver can plan with m -> 0, making T/m blow up
    % at the next linearisation. Hard lower bound at dry mass is safe.
    for k = 1:N
        C = [C; x(7,k) >= m_dry_nd];
    end

    % --- (9) Sigma bounds ---
    C = [C; sigma_low <= sigma <= sigma_high];

    % --- (10) Trust region SOC (per node) ---
    for k = 1:N
        diff_vec = [x(:,k)     - x_prev(:,k);
                    u(:,k)     - u_prev(:,k);
                    Gamma(k)   - Gamma_prev(k)];
        C = [C; norm(diff_vec) <= delta_max + s_delta(k)];
        C = [C; s_delta(k) >= 0];
    end
    % --- (8.5) Monotonic vertical descent (anti-climb) ---

    for k = 2:N
        C = [C; x(6,k) <= 0];
    end
    
    % Thrust regularisation: norm() is SOC-native (YALMIP optimizer can
    % reformulate at compile time), unlike sum(Γ.^2) which is quadratic
    % and rejected by ECOS at runtime. Same minimiser set as the squared
    % version, just a milder penalty rate at the same weight — bump the
    % weight up correspondingly.
    %
    % Effect: prevents bang-coast-bang by penalising the L2 norm of the
    % Γ trajectory. Bang profiles concentrate Γ at the bounds (high norm);
    % smooth profiles spread Γ near the average (lower norm). The optimiser
    % trades a small fuel-suboptimality for a much smoother trajectory.
    %
    % Same convex curvature trick fixes the cone-vertex degeneracy that
    % was causing ECOS exit_code 4 on opening calls.
    w_thrust_reg = 0.02;       % tune in [0.01, 0.2] — too high makes it lazy

    J = -x(7, N) ...
        + w_nu         * norm(nu(:), 1) ...
        + w_delta      * sum(s_delta) ...
        + w_t          * sigma ...
        + w_thrust_reg * norm(Gamma);


    %  Step 6 - Solver options
    options = sdpsettings('solver',  'ecos', ...
                          'verbose', 0,      ...
                          'cachesolvers', 1);


    %  Step 7 - Compile the optimizer
    params_in = {Abar, Bbar, Cbar, Sigbar, zbar, ...
                 x0,   x_target, x_prev,  u_prev, Gamma_prev, ...
                 sigma_low, sigma_high, w_t};

    vars_out  = {x, u, Gamma, sigma, nu, s_delta};

    fprintf('build_socp: compiling SOCP (N=%d, %d intervals)...\n', N, N-1);
    t_compile = tic;
    optimizer_obj = optimizer(C, J, options, params_in, vars_out);
    t_compile = toc(t_compile);
    fprintf('build_socp: compile finished in %.2f s.\n', t_compile);
end