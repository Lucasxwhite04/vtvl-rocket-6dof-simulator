function [theta, phi, psi] = quat2euler(q)
% QUAT2EULER  Convert unit quaternion to ZYX Euler angles.
%
% Used at two interfaces only (Build Plan v2 §1.2, §8.1):
%   1. Inner attitude loop  — consumes theta (pitch) and phi (roll)
%   2. Yaw monitor          — consumes psi (yaw)
% L. White, University of Sheffield, May 2026

% Unpack the 4x1 vector assuming the scalar (real) part is first.
q0 = q(1);   % scalar part
q1 = q(2);   % vector x
q2 = q(3);   % vector y
q3 = q(4);   % vector z

% Calculate the length of the 4D vector and divide each component by it.
q_norm = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
if q_norm > eps
    q0 = q0 / q_norm;
    q1 = q1 / q_norm;
    q2 = q2 / q_norm;
    q3 = q3 / q_norm;
end

% --- ZYX Euler angles from quaternion ---
% Convert the 4D hyper-sphere coordinates into 3D rotational angles.

% Compute Pitch (theta) and clamp the inner math between -1 and 1. 
theta = asin ( max(-1, min(1,  2*(q0*q2 - q3*q1)                        )));

% Compute Roll (phi) and Yaw (psi) using the four-quadrant arctangent.
phi   = atan2(              2*(q0*q1 + q2*q3),   1 - 2*(q1^2 + q2^2)    );
psi   = atan2(              2*(q0*q3 + q1*q2),   1 - 2*(q2^2 + q3^2)    );

end   % quat2euler