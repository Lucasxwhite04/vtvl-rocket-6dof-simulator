function [is_valid, fault_code, fallback_flag, rejection_count_new] = ...
    validity_check(sol, rejection_count, params)
% VALIDITY_CHECK  Pre-splice gate for Mode B solutions.
%
% Examines a candidate SCP solution and decides whether it is safe to
% splice into the active trajectory buffer. Failing solutions are rejected
% and the rejection counter advances. After params.scp.max_rejections
% consecutive rejections (default 3), fallback_flag goes true and the
% guidance subsystem reverts to Mode A for the rest of the phase.
%
% Fail-fast: returns on the first failed check. fault_code identifies
% which one. The order of checks is intentional — cheaper, more reliable
% checks first.
%
%   code 0  - solution is valid
%   code 1  - NaN or Inf anywhere in the state/control trajectory
%   code 2  - LCvx slack exceeds tolerance (convexification non-tight)
%   code 3  - Terminal box: position or velocity error too large
%   code 4  - SCP convergence: nu_norm too large after i_max iterations
%   code 5  - Solver exit code non-zero
%
% --- Inputs ---
%   sol               candidate solution struct, with at least these fields:
%                       .x_k          7 x N   non-dim state trajectory
%                       .u_k          3 x N   non-dim thrust vectors
%                       .Gamma_k      1 x N   non-dim LCvx slack
%                       .sigma        scalar  non-dim time dilation
%                       .x_target     6 x 1   non-dim [r_target; v_target]
%                       .lcvx_slack_N scalar  max LCvx slack (N), SI
%                       .nu_norm      scalar  inf-norm of nu, non-dim
%                       .exit_code    scalar  YALMIP exit code (0=opt)
%   rejection_count   integer, consecutive rejections so far (this phase)
%   params            vehicle params struct
%
% --- Outputs ---
%   is_valid              true if all checks pass
%   fault_code            see above
%   fallback_flag         true if the (potentially incremented)
%                         rejection_count reaches max_rejections
%   rejection_count_new   updated rejection counter
%
% CALLED BY: mpc_solve.m
% CALLS:     nothing
% L. White, University of Sheffield, May 2026

    % --- Check 1: NaN / Inf anywhere in the trajectory ---
    if any(~isfinite(sol.x_k(:))) || any(~isfinite(sol.u_k(:))) ...
            || any(~isfinite(sol.Gamma_k(:))) || ~isfinite(sol.sigma)
        is_valid   = false;
        fault_code = 1;
        rejection_count_new = rejection_count + 1;
        fallback_flag = (rejection_count_new >= params.scp.max_rejections);
        return;
    end

    % --- Check 2: LCvx losslessness ---
    % If slack is large, the relaxation isn't tight - the thrust vector
    % doesn't have magnitude equal to Gamma, meaning the planner is
    % paying for thrust it isn't using. Solution is suspect.
    if sol.lcvx_slack_N > params.scp.lcvx_tol
        is_valid   = false;
        fault_code = 2;
        rejection_count_new = rejection_count + 1;
        fallback_flag = (rejection_count_new >= params.scp.max_rejections);
        return;
    end

    % --- Check 3: Terminal box (position + velocity) ---
    % The SOCP enforces x_k(1:6, N) == x_target as an equality, so this
    % check is normally redundant, any numerical drift would be at the
    % 1e-12 level. But if the solver reported a degraded answer, this
    % catches it.
    r_target_nd = sol.x_target(1:3);
    v_target_nd = sol.x_target(4:6);
    r_end_nd = sol.x_k(1:3, end);
    v_end_nd = sol.x_k(4:6, end);

    pos_err_m   = norm(r_end_nd - r_target_nd) * params.U_L;     % to metres
    vel_err_mps = norm(v_end_nd - v_target_nd) * params.U_V;     % to m/s

    if pos_err_m > params.terminal_pos_tol || ...
       vel_err_mps > params.terminal_vel_tol
        is_valid   = false;
        fault_code = 3;
        rejection_count_new = rejection_count + 1;
        fallback_flag = (rejection_count_new >= params.scp.max_rejections);
        return;
    end

    % --- Check 4: SCP convergence (virtual control) ---
    % Allow some headroom over nu_tol since real cases may not reach
    % machine precision. 10x is the build-plan recommendation.
    if sol.nu_norm > 10 * params.scp.nu_tol
        is_valid   = false;
        fault_code = 4;
        rejection_count_new = rejection_count + 1;
        fallback_flag = (rejection_count_new >= params.scp.max_rejections);
        return;
    end

    % --- Check 5: Solver exit code ---
    % 0 = optimal (YALMIP). Non-zero = something went wrong; even if
    % the trajectory looks fine, we don't trust a solver that flagged
    % an issue.
    if sol.exit_code ~= 0
        is_valid   = false;
        fault_code = 5;
        rejection_count_new = rejection_count + 1;
        fallback_flag = (rejection_count_new >= params.scp.max_rejections);
        return;
    end

    % --- All checks passed: accept and reset rejection counter ---
    is_valid            = true;
    fault_code          = 0;
    rejection_count_new = 0;            % consecutive count resets on accept
    fallback_flag       = false;
end