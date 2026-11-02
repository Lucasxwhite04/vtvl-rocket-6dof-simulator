function [Abar, Bbar, Cbar, Sigbar, zbar] = ...
    linearise_dynamics(x_bar, u_bar, Gamma_bar, sigma_bar, params)
% LINEARISE_DYNAMICS  Discretise the LTV translational dynamics for SCP.
%
% Computes the FOH-discretised Jacobians of the sigma-scaled, non-dimensional
% 3-DoF translational dynamics about a reference iterate, by integrating the
% State Transition Matrix (STM) over each of the N-1 intervals. This is the
% inner machinery the SOCP uses as its dynamics equality constraints.
%
% After Lossless Convexification the control is
% augmented to u_aug = [T_vec; Gamma] in R^4, where Gamma is the slack
% variable satisfying ||T_vec|| <= Gamma. This file produces a SINGLE 7x4
% Jacobian Bbar (concatenating thrust and slack columns), which keeps the
% SOCP construction in build_socp.m simple.
%
% --- Continuous dynamics being linearised ---
%   dx/dtau = sigma * f(x, T_vec, Gamma),     tau in [0, 1] (sigma-scaled)
%
%   f = [ v                       ;
%         g_nd + T_vec/m + D(v)/m ;
%         -Gamma / c_nd            ];          (non-dim)
%
%   x   = [r; v; m] in R^7
%   u   = [T_vec; Gamma] in R^4
%
% --- Linearisation about (x_bar, u_bar, Gamma_bar, sigma_bar) ---
%   A   = sigma_bar * df/dx           7x7    (drag and mass-ratio terms)
%   B   = sigma_bar * df/du_aug       7x4    (I/m in vel row, -1/c in mass row)
%   S   = f(xbar, ubar, Gammabar)     7x1    (df/d(sigma*1) = f)
%   z   = -A * x_bar - B * u_aug_bar  7x1    (affine remainder)
%
% Sanity check on z: sigma*f appears in two places that cancel exactly.
% At the linearisation point we must recover dx/dtau = sigma_bar * fbar:
%   A*xbar + B*ubar + S*sigma_bar + z
%     = A*xbar + B*ubar + fbar*sigma_bar + (-A*xbar - B*ubar)
%     = fbar * sigma_bar    OK
%
% --- FOH discretisation over interval [tau_k, tau_{k+1}] ---
%   x_{k+1} = Abar_k*x_k + Bbar_k*u_k + Cbar_k*u_{k+1}
%             + Sigbar_k*sigma + zbar_k + nu_k
%
%   u(tau) = lam_k(tau)*u_k + lam_{k+1}(tau)*u_{k+1}     (linear interp)
%     lam_k     = (tau_{k+1} - tau) / dtau
%     lam_{k+1} = (tau - tau_k)     / dtau
%
% By Duhamel's principle, with Phi(tau, xi) the STM:
%   Abar_k  = Phi(tau_{k+1}, tau_k)
%   Bbar_k  = integral of Phi(tau_{k+1}, xi) * B(xi) * lam_k(xi)     dxi
%   Cbar_k  = integral of Phi(tau_{k+1}, xi) * B(xi) * lam_{k+1}(xi) dxi
%   Sigbar  = integral of Phi(tau_{k+1}, xi) * S(xi)                 dxi
%   zbar    = integral of Phi(tau_{k+1}, xi) * z(xi)                 dxi
%
% Computed WITHOUT matrix inverses by integrating the augmented ODE:
%   d/dtau [Phi;      = [ A*Phi;            Phi(tau_k)   = I_7
%           Bint;       A*Bint + lam_k*B;   Bint(tau_k)  = 0
%           Cint;       A*Cint + lam_{k+1}*B; Cint(...)  = 0
%           Sint;       A*Sint + S;         Sint(tau_k)  = 0
%           zint]       A*zint + z]         zint(tau_k)  = 0
%
% This follows directly from variation of constants applied to each
% column of the convolution integral.
%
% --- Dispersions policy ---
% This function uses nominal drag and exhaust-velocity values. It does
% NOT apply params.dispersions.scale.* — those represent the truth that
% the MPC cannot know. The SCP virtual control nu_k absorbs the
% nominal-vs-truth gap during closed-loop operation.
%
% --- Inputs ---
%   x_bar      - reference state iterate       [non-dim]   7 x N
%   u_bar      - reference thrust vector       [non-dim]   3 x N
%   Gamma_bar  - reference LCvx slack          [non-dim]   1 x N
%   sigma_bar  - reference time dilation       [non-dim]   scalar
%   params     - vehicle params struct
%
% --- Outputs (each a cell array of length N-1) ---
%   Abar{k}    - discretised state Jacobian        7 x 7
%   Bbar{k}    - discretised augmented-control Jacobian at start   7 x 4
%                 columns 1:3 -> T_vec   column 4 -> Gamma
%   Cbar{k}    - discretised augmented-control Jacobian at end     7 x 4
%   Sigbar{k}  - discretised sigma column                          7 x 1
%   zbar{k}    - discretised affine remainder                      7 x 1
%
% CALLED BY: scp_iteration.m
% CALLS:     nothing (self-contained; local function aug_ode below)
% L. White, University of Sheffield, May 2026

    % --- Shape and sanity checks ---
    N = params.scp.N;
    assert(size(x_bar, 1) == 7 && size(x_bar, 2) == N, ...
           'x_bar must be 7 x %d (got %d x %d)', N, size(x_bar,1), size(x_bar,2));
    assert(size(u_bar, 1) == 3 && size(u_bar, 2) == N, ...
           'u_bar must be 3 x %d', N);
    assert(numel(Gamma_bar) == N, 'Gamma_bar must have N=%d entries', N);
    assert(isscalar(sigma_bar), 'sigma_bar must be scalar');

    Gamma_bar = Gamma_bar(:).';  % force 1 x N row

    % --- Non-dim constants ---
    c_nd = params.c / params.U_V;     % non-dim exhaust velocity

    % Drag coefficients folded into single non-dim scalars:
    %   D_nd_i(v_nd) = -k_i * v_nd_i * |v_nd_i|
    
    % Derivation: D_dim = -0.5*rho*CD*A * v_dim*|v_dim|
    %             = -0.5*rho*CD*A * (v_nd*U_V) * |v_nd*U_V|
   
    % In non-dim form D_nd = D_dim / U_F, so coefficient is
    %   k = 0.5 * rho * CD * A * U_V^2 / U_F
    kx = 0.5 * params.rho * params.CDxy * params.A_lateral * ...
         params.U_V^2 / params.U_F;
    kz = 0.5 * params.rho * params.CDz  * params.A_axial   * ...
         params.U_V^2 / params.U_F;

    % Pack into struct for the local ODE function (avoids long arg list)
    L.c_nd = c_nd;
    L.kx   = kx;
    L.kz   = kz;
    L.sigma = sigma_bar;

    % --- Discretisation grid ---
    dtau = 1.0 / (N - 1);

    % --- Augmented ODE: state layout ---
    %   y(1:49)   = vec(Phi)        7x7
    %   y(50:77)  = vec(Bint)       7x4
    %   y(78:105) = vec(Cint)       7x4
    %   y(106:112)= Sint            7x1
    %   y(113:119)= zint            7x1
    n_Phi = 49;  n_B = 28;  n_C = 28;  n_S = 7;  n_z = 7;
    n_total = n_Phi + n_B + n_C + n_S + n_z;     % = 119

    % Initial condition for each interval: Phi = I, all integrals = 0
    y0 = zeros(n_total, 1);
    y0(1:n_Phi) = reshape(eye(7), [], 1);

    % --- Allocate outputs ---
    Abar   = cell(N-1, 1);
    Bbar   = cell(N-1, 1);
    Cbar   = cell(N-1, 1);
    Sigbar = cell(N-1, 1);
    zbar   = cell(N-1, 1);

    % --- ODE options ---
    % Tight tolerances: the linearisation accuracy directly affects SCP
    % convergence. RelTol 1e-8 is well within ODE45's stable range and
    % adds only milliseconds vs the default 1e-3.
    ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

    % --- Loop over intervals ---
    for k = 1:N-1
        tau_k   = (k - 1) * dtau;
        tau_kp1 =  k      * dtau;

        % Endpoint reference values for this interval
        L.x_k    = x_bar(:, k);
        L.x_kp1  = x_bar(:, k+1);
        L.T_k    = u_bar(:, k);
        L.T_kp1  = u_bar(:, k+1);
        L.G_k    = Gamma_bar(k);
        L.G_kp1  = Gamma_bar(k+1);
        L.tau_k  = tau_k;
        L.dtau   = dtau;

        % Integrate
        [~, Y] = ode45(@(tau, y) aug_ode(tau, y, L), ...
                       [tau_k, tau_kp1], y0, ode_opts);
        yf = Y(end, :).';

        % Unpack
        idx = 0;
        Abar{k}   = reshape(yf(idx+1:idx+n_Phi), 7, 7);  idx = idx + n_Phi;
        Bbar{k}   = reshape(yf(idx+1:idx+n_B  ), 7, 4);  idx = idx + n_B;
        Cbar{k}   = reshape(yf(idx+1:idx+n_C  ), 7, 4);  idx = idx + n_C;
        Sigbar{k} =         yf(idx+1:idx+n_S  );         idx = idx + n_S;
        zbar{k}   =         yf(idx+1:idx+n_z  );
    end
end


%%  LOCAL FUNCTION — augmented ODE for STM + convolution integrals
function dy = aug_ode(tau, y, L)
% AUG_ODE  Right-hand side of the augmented STM-integration ODE.
%
% Computes A(tau), B_aug(tau), S(tau), z(tau) along the linearisation path,
% then assembles the derivative of the stacked state vector y.

    % --- Interpolation weights for the reference iterate ---
    lam_k   = (L.tau_k + L.dtau - tau) / L.dtau;   % weight on node k
    lam_kp1 = (tau - L.tau_k)          / L.dtau;   % weight on node k+1

    % Reference values at this tau (linear interp inside interval)
    x   = lam_k * L.x_k   + lam_kp1 * L.x_kp1;
    Tvec= lam_k * L.T_k   + lam_kp1 * L.T_kp1;
    G   = lam_k * L.G_k   + lam_kp1 * L.G_kp1;

    v = x(4:6);
    m = max(x(7), 0.05);    % guard against pathological zero-mass

    % --- Drag and its Jacobian (non-dim) ---
    %   D_i = -k_i * v_i * |v_i|     so   d/dv_i (v_i*|v_i|) = 2|v_i|
    D = [-L.kx * v(1) * abs(v(1));
         -L.kx * v(2) * abs(v(2));
         -L.kz * v(3) * abs(v(3))];

    dD_dv = diag([-L.kx * 2*abs(v(1));
                  -L.kx * 2*abs(v(2));
                  -L.kz * 2*abs(v(3))]);    % 3x3 diagonal

    % --- Build the Jacobian matrices in PHYSICAL (un-scaled) form ---
    % A_phys = df/dx (7x7), B_phys = df/du_aug (7x4), then multiply by sigma.
    g_nd = [0; 0; -1];                       % non-dim gravity

    % Non-dim dynamics value f
    f = [v;
         g_nd + Tvec/m + D/m;
         -G / L.c_nd];

    % df/dx blocks
    A_phys = zeros(7, 7);
    A_phys(1:3, 4:6) = eye(3);                       % dr/dot/dv
    A_phys(4:6, 4:6) = dD_dv / m;                    % dvdot/dv (drag)
    A_phys(4:6, 7)   = -(Tvec + D) / m^2;            % dvdot/dm
    % Mass row is zero (d/dx of -G/c is zero w.r.t. state)

    % df/du_aug = [df/dT_vec | df/dGamma]  (7x4)
    B_phys = zeros(7, 4);
    B_phys(4:6, 1:3) = eye(3) / m;                   % dvdot/dT_vec
    B_phys(7,   4)   = -1 / L.c_nd;                  % dmdot/dGamma

    % Scale by sigma to form A and B
    A     = L.sigma * A_phys;
    B_aug = L.sigma * B_phys;

    % df/d(sigma*1) = f at reference (NOT multiplied by sigma)
    S = f;

    % Affine remainder z (see header for cancellation derivation)
    u_aug = [Tvec; G];
    z = -A*x - B_aug*u_aug;

    % --- Unpack the augmented state ---
    n_Phi = 49;  n_B = 28;  n_C = 28;  n_S = 7;  n_z = 7;
    idx = 0;
    Phi   = reshape(y(idx+1:idx+n_Phi), 7, 7);  idx = idx + n_Phi;
    Bint  = reshape(y(idx+1:idx+n_B  ), 7, 4);  idx = idx + n_B;
    Cint  = reshape(y(idx+1:idx+n_C  ), 7, 4);  idx = idx + n_C;
    Sint  =         y(idx+1:idx+n_S  );         idx = idx + n_S;
    zint  =         y(idx+1:idx+n_z  );

    % --- Augmented derivatives ---
    dPhi  = A * Phi;
    dBint = A * Bint + lam_k   * B_aug;
    dCint = A * Cint + lam_kp1 * B_aug;
    dSint = A * Sint + S;
    dzint = A * zint + z;

    % --- Repack as column vector ---
    dy = [reshape(dPhi,  [], 1);
          reshape(dBint, [], 1);
          reshape(dCint, [], 1);
          dSint;
          dzint];
end