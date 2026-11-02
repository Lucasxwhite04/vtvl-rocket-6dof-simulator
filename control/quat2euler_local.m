function [theta_hat, phi_hat, psi_hat] = quat2euler_local(q)
% QUAT2EULER_LOCAL  Convert unit quaternion to ZYX Euler angles (pitch, roll, yaw).
%
% Converts the MEKF quaternion state (indices 7:10 of the 17-state vector)
% to Euler angles for downstream consumers.
% L. White May 2026

% Unpack the quaternion assuming the scalar (real) part is first.
qw = q(1); qx = q(2); qy = q(3); qz = q(4);

%%  DEFENSIVE RENORMALIZATION
% Calculate the length of the 4D vector and divide each component by it.
q_norm = sqrt(qw^2 + qx^2 + qy^2 + qz^2);
if q_norm > eps
    qw = qw / q_norm;
    qx = qx / q_norm;
    qy = qy / q_norm;
    qz = qz / q_norm;
end

%%  EULER ANGLE EXTRACTION (ZYX Convention)
% --- Pitch (theta) ---
% Compute Pitch and clamp the inner math strictly between -1 and 1.
theta_hat = asin( max(-1, min(1,  2*(qw*qy - qz*qx) )) );

% --- Roll (phi) ---
% Compute Roll using the four-quadrant arctangent (atan2).
phi_hat = atan2( 2*(qw*qx + qy*qz), ...
                 1 - 2*(qx^2 + qy^2) );

% --- Yaw (psi) ---
% Compute Yaw using the four-quadrant arctangent.
psi_hat = atan2( 2*(qw*qz + qx*qy), ...
                 1 - 2*(qy^2 + qz^2) );

end