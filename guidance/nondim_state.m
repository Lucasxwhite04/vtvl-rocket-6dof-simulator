function y = nondim_state(x, params, direction)
% NONDIM_STATE  Convert a 7-state [r; v; m] between SI and non-dimensional units.
%
% Mode B's SCP solves in non-dimensional units to keep the SOCP well
% conditioned. Every file in guidance/ that touches
% state vectors calls this function exactly once on entry (to non-dim) and
% exactly once on exit (to re-dim), so that the internal arithmetic always
% works in the conditioned domain.
%
% Scaling factors (defined in plant_params.m):
%   U_L = 100 m          length scale
%   U_T ≈ 3.19 s         time scale,   sqrt(U_L / g0)
%   U_V ≈ 31.3 m/s       velocity = U_L / U_T
%   U_M = 16.88 kg       mass scale   = m0
%
% INPUTS:
%   x         - 7-state [r(3); v(3); m] in current units             7x1 column
%               OR (7 x N) matrix of N states stored column-wise
%   params    - vehicle params struct (must contain U_L, U_V, U_M)
%   direction - 'to_nd'  : SI    → non-dim
%             - 'to_si'  : non-dim → SI
%
% OUTPUTS:
%   y         - same size as x, in the requested units
%
% CALLED BY: mpc_solve.m, scp_iteration.m, warm_start.m, extract_command.m
% CALLS:     nothing
% L. White, University of Sheffield, May 2026

% --- Defensive size check ---
if size(x, 1) ~= 7
    error('nondim_state:badSize', ...
        'Input must have 7 rows (got %d).', size(x, 1));
end

% --- Scaling vector - applied row-wise ---
% Rows 1–3: position    - divide by U_L
% Rows 4–6: velocity    - divide by U_V
% Row    7: mass        - divide by U_M
scale = [params.U_L; params.U_L; params.U_L;
    params.U_V; params.U_V; params.U_V;
    params.U_M];   % (7 x 1)

switch lower(direction)
    case 'to_nd'
        y = x ./ scale;        % implicit expansion across columns
    case 'to_si'
        y = x .*  scale;
    otherwise
        error('nondim_state:badDir', ...
            'direction must be ''to_nd'' or ''to_si'', got "%s".', ...
            direction);
end
end