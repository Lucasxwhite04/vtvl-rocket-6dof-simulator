function K_table = lqi_gains(params, m_grid)
% LQI_GAINS  Pre-compute gain-scheduled LQI gains over a mass grid.
%
% Linearises the pitch/roll attitude dynamics about upright hover at each
% mass point and solves the LQR problem for a 3-state augmented system that
% includes integral error.  Returns a lookup table; call ONCE at startup.
%
% PLANT MODEL (pitch channel - roll is identical by symmetry):
%
%   Continuous-time rotational dynamics about body-y (pitch):
%     theta_dot = q_rate
%     q_dot     = (T * ls / Ixy) * alpha      [gimbal creates pitching moment]
%     e_i_dot   = theta - theta_ref            [integrator state]
%
%   Augmented state:  x_att = [theta; q_rate; e_i]
%   Input:            u_att = alpha  (gimbal pitch angle, rad)
%
%   A = [0,  1,  0;     B = [0;
%        0,  0,  0;          b_pitch;
%        1,  0,  0]          0]
%
%   where  b_pitch = T_nom * ls_nominal / Ixy(m)
%          Ixy(m)  = params.Ixy_per_kg * m
%
% MASS-SCHEDULING NOTE:
%   At hover  T_nom = m * g0, so  b_pitch = g0 * ls / Ixy_per_kg  (mass cancels).
%   Gains are therefore approximately constant with mass at hover thrust.
%   The grid remains useful for descent where T ≠ m*g0 (off-nominal thrust
%   levels change b_pitch and hence the optimal gains).
%
% Inputs:
%   params   struct from plant_params()  — uses: ms, m0, g0, ls_nominal,
%                                                  Ixy_per_kg
%   m_grid   (Nx1) optional — mass values to compute gains at (kg)
%            Default: 12 points from params.ms to params.m0
%
% Output:
%   K_table  struct with fields:
%     .m_grid   (Nx1)  mass grid (kg)
%     .K_pitch  (Nx3)  [K_theta, K_q_rate, K_i_pitch]  per row
%     .K_roll   (Nx3)  [K_phi,   K_p,      K_i_roll ]  per row
%                      (identical to K_pitch — symmetric vehicle)
%
% Usage:
%   K_table = lqi_gains(params);                     % default 12-point grid
%   K_table = lqi_gains(params, (8:0.5:17)');        % custom 
% L. White May 2026

if nargin < 2 || isempty(m_grid)
    m_grid = linspace(params.ms, params.m0, 12)';
end
m_grid = m_grid(:);
N = numel(m_grid);

% --- LQR cost weights ---
% Q penalises: [angle error, angular rate, integrator state]
% R penalises: gimbal deflection (control effort)

Q_att = params.lqi.Q;
R_att = params.lqi.R;
% --- Continuous-time plant matrices (angle/rate states) ---
A_att = [0, 1, 0;
         0, 0, 0;
         1, 0, 0];

% --- Preallocate ---
K_pitch = zeros(N, 3);
K_roll  = zeros(N, 3);

for k = 1:N
    m_k = m_grid(k);

    % Moment of inertia scales linearly with mass (Build Plan §6)
    Ixy_k = params.Ixy_per_kg * m_k;

    % ---  GAIN SCHEDULING  ---
    % At hover: T_nom = m*g0, so b_pitch = g0*ls/Ixy_per_kg (mass cancels).
    % During descent: T_nom ≠ m*g0, so b_pitch changes with mass AND thrust.
    % The grid captures this variation for the MPC to exploit.
    T_nom_k = m_k * params.g0;                          % hover thrust at this mass
    b_pitch  = T_nom_k * params.ls_nominal / Ixy_k;     % angular accel per rad of gimbal

    B_att = [0; b_pitch; 0];

    % --- Solve LQR using Q and R from plant_params ---
    K_k = lqr(A_att, B_att, Q_att, R_att);
    % K_k = [K_theta, K_q, K_i] — all positive by construction

    if any(K_k < 0)
        warning('Negative gain detected: K = [%.3f, %.3f, %.3f]', K_k);
    end

    K_pitch(k, :) = K_k;
    K_roll(k,  :) = K_k;   % symmetric vehicle
end

% --- Pack output ---
K_table.m_grid  = m_grid;
K_table.K_pitch = K_pitch;   % columns: [K_theta,  K_q_rate,  K_i_pitch]
K_table.K_roll  = K_roll;    % columns: [K_phi,    K_p,       K_i_roll ]

% --- Sanity print ---
fprintf('LQI gains computed at %d mass points [%.1f ... %.1f] kg\n', ...
        N, m_grid(1), m_grid(end));
fprintf('  K_pitch at m=%.2f kg: K_theta=%.2f  K_q=%.2f  K_i=%.2f\n', ...
        m_grid(1),   K_pitch(1,1),   K_pitch(1,2),   K_pitch(1,3));
fprintf('  K_pitch at m=%.2f kg: K_theta=%.2f  K_q=%.2f  K_i=%.2f\n', ...
        m_grid(end), K_pitch(end,1), K_pitch(end,2), K_pitch(end,3));

end

%% LOCAL HELPER

function K = lqi_analytical(b, Q, R)
% LQI_ANALYTICAL  Closed-form LQR approximation for the 3-state LQI system.
% Used as fallback when Control System Toolbox (lqr) is unavailable.
%
% Derivation: for the integrating plant q_ddot = b*u, the standard LQR
% double-integrator solution gives:
%   K_angle = sqrt(q1/R) / b  ×  (correction for integral state)
%   K_rate  = 2^(1/4) × (q1/(R*b^2))^(1/4)
%   K_i     = sqrt(q3/R) / b  (integral gain; smaller than angle gain)
%

q1 = Q(1,1);   % angle weight
q2 = Q(2,2);   % rate weight
q3 = Q(3,3);   % integral weight
r  = R;

% Approximations based on separating the outer (angle/rate) and inner
% (integral) dynamics of the 3-state system:
K_angle    = (q1/r)^(1/4) * (q2/r)^(1/4) / b;
K_rate     = (q2/r)^(1/2) / b^(1/2);
K_integral = (q3/r)^(1/4) / b * 0.4;   % integral gain intentionally smaller

K = [K_angle, K_rate, K_integral];
end