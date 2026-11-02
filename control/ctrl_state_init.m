function s = ctrl_state_init(m0, theta0, phi0)
% CTRL_STATE_INIT  Initialise attitude controller persistent state.
%
% theta0 and phi0 default to zero (upright launch condition).
% Call once at startup or on phase transition reset; pass the result
% as ctrl_in on the first call to attitude_control().
%
% Inputs:
%   m0      (scalar)  initial vehicle mass (kg)
%   theta0  (rad)     initial pitch angle (default 0)
%   phi0    (rad)     initial roll  angle (default 0)
%
% Output:
%   s   struct with all persistent controller fields.
%
% USAGE:
%   ctrl = ctrl_state_init(params.m0);           % upright launch
%   ctrl = ctrl_state_init(params.m0, 0, 0);     % explicit upright
%   ctrl = ctrl_state_init(params.m0, th, ph);   % non-zero initial attitude
% L. White May 2026

if nargin < 2, theta0 = 0; end
if nargin < 3, phi0   = 0; end

s.i_theta      = 0;       % pitch integrator state  (rad·s)
s.i_phi        = 0;       % roll  integrator state  (rad·s)
s.theta_prev   = theta0;  % previous pitch angle    (rad) - seeds rate estimator
s.phi_prev     = phi0;    % previous roll  angle    (rad) - seeds rate estimator
s.e_theta_prev = 0;       % previous pitch error    (rad) - for diagnostics
s.e_phi_prev   = 0;       % previous roll  error    (rad) - for diagnostics
s.m_est        = m0;      % initial mass            (kg)  - for diagnostics
end