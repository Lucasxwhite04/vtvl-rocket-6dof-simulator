function [x_new, u_new, Gamma_new, sigma_new, info] = scp_iteration( ...
    x_bar, u_bar, Gamma_bar, sigma_bar, ...
    x0, x_target, sigma_low, sigma_high, w_t, ...
    optimizer_obj, params)
% SCP_ITERATION  One outer-loop iteration of Successive Convex Programming.
%
% Step sequence:
%   1. Linearise the 3-DoF dynamics about (x_bar, u_bar, Gamma_bar, sigma_bar)
%      to produce the per-interval discretised Jacobians.
%   2. Pack them into the 3D-array form the compiled SOCP expects.
%   3. Call the compiled YALMIP optimizer with x_prev = x_bar etc.
%   4. Unpack the solution and compute convergence metrics:
%        - nu_norm    : inf-norm of virtual control across all intervals
%        - delta_norm : inf-norm of change in (x, u, Gamma) from iterate k
%        - lcvx_slack : max |Gamma - ||u_T||| in Newtons
%
% --- Inputs (all non-dim unless noted) ---
%   x_bar       - state iterate          7 x N
%   u_bar       - thrust vector iterate  3 x N
%   Gamma_bar   - LCvx slack iterate     1 x N
%   sigma_bar   - time dilation iterate  scalar
%   x0          - initial condition      7 x 1
%   x_target    - target [r; v]          6 x 1
%   sigma_low   - sigma lower bound      scalar  (= sigma_high for landing)
%   sigma_high  - sigma upper bound      scalar
%   w_t         - time penalty weight    scalar  (0 landing, 1 descent)
%   optimizer_obj - compiled YALMIP optimizer from build_socp.m
%   params      - vehicle params struct
%
% --- Outputs ---
%   x_new, u_new, Gamma_new, sigma_new - new iterate (non-dim)
%   info       - struct with the following fields:
%                  .exit_code   YALMIP solver exit code (0 = optimal)
%                  .nu_norm     max ||nu_k|| over intervals      [non-dim]
%                  .delta_norm  max change in iterate (inf-norm) [non-dim]
%                  .lcvx_slack  max |Gamma - ||u_T|||            [N, SI]
%                  .s_delta_max max trust-region slack           [non-dim]
%                  .nu          full virtual-control matrix      7 x (N-1)
%                  .s_delta     full trust-region slack vector   1 x N
%                  .t_lin       linearisation wall-clock         [s]
%                  .t_solve     SOCP solve wall-clock            [s]
%
% CALLED BY: mpc_solve.m, test_modeB_step4.m
% CALLS:     linearise_dynamics.m, the compiled optimizer
%
% L. White, University of Sheffield, May 2026

    %% Step 1 - Linearise about current iterate
    t0 = tic;
    [Abar_c, Bbar_c, Cbar_c, Sigbar_c, zbar_c] = ...
        linearise_dynamics(x_bar, u_bar, Gamma_bar, sigma_bar, params);
    t_lin = toc(t0);

    %% Step 2 - Repack cells into 3D arrays for the optimizer
    Abar3D   = cat(3, Abar_c{:});      % 7 x 7 x (N-1)
    Bbar3D   = cat(3, Bbar_c{:});      % 7 x 4 x (N-1)
    Cbar3D   = cat(3, Cbar_c{:});      % 7 x 4 x (N-1)
    Sigbar2D = [Sigbar_c{:}];          % 7 x (N-1)
    zbar2D   = [zbar_c{:}];            % 7 x (N-1)

    %% Step 3 - Call the compiled optimizer
    % Parameter order MUST match build_socp.m
    param_values = {Abar3D, Bbar3D, Cbar3D, Sigbar2D, zbar2D, ...
                    x0, x_target, x_bar, u_bar, Gamma_bar, ...
                    sigma_low, sigma_high, w_t};

    t0 = tic;
    [sol, exit_code] = optimizer_obj(param_values);
    t_solve = toc(t0);

    if iscell(exit_code), exit_code = exit_code{1}; end

    %% Step 4 - Unpack solution
    % If the solver failed, return the previous iterate untouched and
    % flag a large nu_norm to force the outer loop to bail / fall back
    if exit_code ~= 0
        x_new     = x_bar;
        u_new     = u_bar;
        Gamma_new = Gamma_bar;
        sigma_new = sigma_bar;
        info = struct('exit_code', exit_code, ...
                      'nu_norm',    inf, ...
                      'delta_norm', inf, ...
                      'lcvx_slack', inf, ...
                      's_delta_max', inf, ...
                      'nu',         [], ...
                      's_delta',    [], ...
                      't_lin',      t_lin, ...
                      't_solve',    t_solve);
        return;
    end

    x_new       = sol{1};       % 7 x N
    u_new       = sol{2};       % 3 x N
    Gamma_new   = sol{3};       % 1 x N
    sigma_new   = sol{4};       % scalar
    nu_sol      = sol{5};       % 7 x (N-1)
    s_delta_sol = sol{6};       % 1 x N

    %% Step 5 - Convergence metrics
    % nu_norm: largest dynamics-defect among intervals. Should -> 0 as the
    % iterate stops moving and the linearisation becomes self-consistent.
    nu_norm = max(vecnorm(nu_sol));       % scalar, inf-norm across k

    % delta_norm: largest per-node change in (x, u, Gamma) from iterate.
    % Use inf-norm of the stacked difference vector
    delta_per_node = zeros(1, params.scp.N);
    for k = 1:params.scp.N
        d = [x_new(:,k)   - x_bar(:,k);
             u_new(:,k)   - u_bar(:,k);
             Gamma_new(k) - Gamma_bar(k)];
        delta_per_node(k) = norm(d, inf);
    end
    delta_norm = max(delta_per_node);

    % LCvx slack in SI units (Newtons) 
    lcvx_slack_nd = max(abs(Gamma_new - vecnorm(u_new)));
    lcvx_slack    = lcvx_slack_nd * params.U_F;     % convert to N

    s_delta_max = max(s_delta_sol);

    %% Step 6 - Pack info
    info = struct( ...
        'exit_code',   exit_code, ...
        'nu_norm',     nu_norm, ...
        'delta_norm',  delta_norm, ...
        'lcvx_slack',  lcvx_slack, ...
        's_delta_max', s_delta_max, ...
        'nu',          nu_sol, ...
        's_delta',     s_delta_sol, ...
        't_lin',       t_lin, ...
        't_solve',     t_solve);
end