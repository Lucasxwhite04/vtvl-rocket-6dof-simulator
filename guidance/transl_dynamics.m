function xdot = transl_dynamics(x, T_vec, Gamma, params)
% TRANSL_DYNAMICS  3-DoF translational dynamics for the MPC and state predictor.
%
% Continuous-time, dimensional (SI). This is the model the online MPC
% plans about. It is not the 6-DoF truth plant, attitude is handled by the 
% inner LQI loop and presumed to have already settled to align body-z with 
% T_vec by the time the next MPC call runs.
%
% After Lossless Convexification (LCvx), the thrust vector T_vec
% and the slack scalar Γ are treated as SEPARATE decision variables:
%   • T_vec    drives the velocity equation (linear in T_vec)
%   • Γ        drives the mass-flow equation (linear in Γ)
%   • coupling: enforced by SOC constraint ||T_vec|| ≤ Γ in build_socp.m
%
% This function evaluates f(x, T_vec, Γ) directly — it does NOT enforce
% the LCvx coupling itself. The caller (state_predictor, linearise_dynamics)
% supplies T_vec and Γ values that already satisfy ||T_vec|| ≈ Γ.
%
% --- State vector (7×1) ---
%   x(1:3)  inertial position r          [m]
%   x(4:6)  inertial velocity v          [m/s]
%   x(7)    total mass m                 [kg]
%
% --- Inputs ---
%   x       - 7×1 state vector (SI)
%   T_vec   - 3×1 thrust vector, inertial frame [N]
%   Gamma   - scalar LCvx slack ≥ ||T_vec|| [N]
%   params  - vehicle params struct (uses: c, grav, rho, CDxy, CDz,
%             A_lateral, A_axial, optionally dispersions.scale)
%
% --- Output ---
%   xdot    - 7×1 state derivative (SI)
%
% --- Drag model ---
% Same smooth C¹ form as plant_dynamics.m:
%   D_i = -0.5 · ρ · CD_i · A_i · v_i · |v_i|
%
% IMPORTANT SIMPLIFICATION: the truth plant resolves drag in body-frame
% components (CDxy on body x/y, CDz on body z), because the body
% orientation determines projected area. The MPC has no attitude state.
% Approximate by using inertial-frame components with an "effective"
% drag, assuming the body is roughly vertical, so:
%   inertial x, y → lateral drag (CDxy, A_lateral)
%   inertial z    → axial   drag (CDz,  A_axial)
% Linearisation error is small for near-vertical descent and is absorbed
% by the SCP virtual-control term ν.
%
% CALLED BY: state_predictor.m, linearise_dynamics.m
% CALLS:     nothing
%
% L. White, University of Sheffield, May 2026
% --- Unpack state ---
v = x(4:6);
m = x(7);

% --- Apply Monte-Carlo dispersion scale factors (nominal = 1) ---
% Same convention as plant_dynamics.m — keeps the predictor's model
% consistent with the truth plant in dispersion runs.
rho  = params.rho  * params.dispersions.scale.rho_scale;
CDxy = params.CDxy * params.dispersions.scale.CDxy_scale;
CDz  = params.CDz  * params.dispersions.scale.CDz_scale;
c    = params.c    * params.dispersions.scale.c_scale;

% --- Drag (inertial-frame approximation, see header note) ---
D = [-0.5 * rho * CDxy * params.A_lateral * v(1) * abs(v(1));
    -0.5 * rho * CDxy * params.A_lateral * v(2) * abs(v(2));
    -0.5 * rho * CDz  * params.A_axial   * v(3) * abs(v(3))];   % 3x1 [N]

% --- Translational kinematics + dynamics ---
%   ṙ = v
%   v̇ = g + T_vec/m + D/m
%   ṁ = -Γ / c          (uses LCvx slack, not ||T_vec||)
xdot = zeros(7, 1);
xdot(1:3) = v;
xdot(4:6) = params.grav + T_vec / m + D / m;
xdot(7)   = -Gamma / c;
end