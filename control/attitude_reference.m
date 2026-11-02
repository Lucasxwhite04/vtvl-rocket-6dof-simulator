function [theta_ref, phi_ref] = attitude_reference(n_hat_cmd)
% ATTITUDE_REFERENCE  Convert inertial thrust direction to body attitude references.
%
% Answers: "if thrust must point in direction n_hat_cmd in the inertial
% frame, what pitch (theta) and roll (phi) must the rocket body be at?"
%
% INPUTS:
%   n_hat_cmd - commanded thrust unit direction  [inertial frame]  3x1
%
% OUTPUTS:
%   theta_ref - commanded pitch angle  [rad]  scalar
%   phi_ref   - commanded roll  angle  [rad]  scalar
%
% CALLED BY: attitude_control.m
% CALLS:     nothing
% L. White May 2026

% Normalise input defensively — guidance should already output unit vector
% but numerical drift can accumulate over long runs
n = n_hat_cmd(:) / norm(n_hat_cmd);

% Recover pitch and roll from the thrust direction vector.
% For an upright rocket, column 3 of R_bi (body Z in inertial frame) is:
%   [sin(theta); -sin(phi)*cos(theta); cos(theta)*cos(phi)]
% Since n_hat_cmd must equal this column:
%   n(1) =  sin(theta)  →  theta = asin( n(1))
%   n(2) = -sin(phi)    →  phi   = asin(-n(2))
%
% Clamp argument to [-1,1] to guard against numerical noise pushing
% slightly outside the valid domain of asin (would return NaN otherwise)
theta_ref = asin(max(-1, min(1,  n(1))));   % pitch [rad]
phi_ref   = asin(max(-1, min(1, -n(2))));   % roll  [rad]

% Body tilt limit - this is NOT the gimbal limit (±5 deg).
% The rocket body can legitimately tilt up to ±30 deg during
% aggressive lateral manoeuvres in ascent and descent.
MAX_ATT   = 30 * pi/180;                    % ±35 deg [rad]
theta_ref = max(-MAX_ATT, min(MAX_ATT, theta_ref));
phi_ref   = max(-MAX_ATT, min(MAX_ATT, phi_ref));

end