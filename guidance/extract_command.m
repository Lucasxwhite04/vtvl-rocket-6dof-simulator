function [T_cmd, n_hat_cmd, lcvx_slack_N] = ...
    extract_command(u_first_nd, Gamma_first_nd, params)
% EXTRACT_COMMAND  Pull the first-step thrust command from the SCP solution.
%
% Reads the thrust vector and LCvx slack at the FIRST node of the converged
% SCP solution, applies the losslessness check, and re-dimensionalises into
% SI units suitable for handoff to attitude_reference.m.
%
% Only the FIRST node is used as the actual command, the rest of the trajectory 
% is stored in the trajectory buffer for the inner attitude loop to read at its 
% own rate. The full solution is what gets spliced; this function only 
% produces the next outgoing command.
%
% --- Inputs (non-dim) ---
%   u_first_nd      3 x 1    thrust vector at node 1 (non-dim force)
%   Gamma_first_nd  1 x 1    LCvx slack at node 1 (non-dim force)
%   params          struct   vehicle params (uses U_F, Tmin, Tmax)
%
% --- Outputs (SI) ---
%   T_cmd           scalar   thrust magnitude command [N]
%                            clamped to [Tmin, Tmax]
%   n_hat_cmd       3 x 1    thrust unit direction (inertial frame)
%                            fallback = [0; 0; 1] if numerical zero
%   lcvx_slack_N    scalar   |Gamma_1 - ||u_1||| in Newtons
%                            for downstream validity_check
%
% CALLED BY: mpc_solve.m
% CALLS:     nothing
%
% L. White, University of Sheffield, May 2026

% --- LCvx losslessness check (non-dim, then convert to N) ---
u_mag_nd     = norm(u_first_nd);
lcvx_slack_nd = abs(Gamma_first_nd - u_mag_nd);
lcvx_slack_N  = lcvx_slack_nd * params.U_F;

% --- Re-dimensionalise thrust magnitude ---
% Use ||u|| rather than Gamma because the actual force on the
% vehicle is determined by the thrust vector itself. If LCvx is
% truly lossless (slack ~ 0), the two are interchangeable; if it
% has degraded, ||u|| is the physically correct one.
T_cmd_raw_N = u_mag_nd * params.U_F;
T_cmd       = max(params.Tmin, min(params.Tmax, T_cmd_raw_N));

% --- Unit thrust direction (frame-invariant, no scaling) ---
% Fallback to body-z up if the thrust vector is numerical noise.
% This is purely defensive - LCvx box forces Gamma >= Tmin > 0, so
% a sane solution should never trigger this branch.
if u_mag_nd > 1e-9
    n_hat_cmd = u_first_nd / u_mag_nd;
else
    n_hat_cmd = [0; 0; 1];
end
end