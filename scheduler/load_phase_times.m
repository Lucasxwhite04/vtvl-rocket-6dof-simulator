function t_phase_nominal = load_phase_times(ref_traj_path, use_hover)
% LOAD_PHASE_TIMES  Extract nominal phase end times from reference trajectory.
%
% Returns a 3-element vector when hover is disabled, or 4-element
% when hover is enabled.
%
%   USE_HOVER = false  -  [t1_end, t3_end, t4_end]       (3x1)
%   USE_HOVER = true   -  [t1_end, t2_end, t3_end, t4_end] (4x1)
%
% Inputs:
%   ref_traj_path  (char)   path to reference_trajectory.mat
%                           Pass '' or omit to use fallback nominal times.
%   use_hover      (bool)   match USE_HOVER flag in phase_scheduler.m
%                           Default: false
%
% Output:
%   t_phase_nominal  (3x1 or 4x1)  nominal phase end times (s)
% L. White, University of Sheffield, May 2026


if nargin < 2
    use_hover = false;
end

% Fallback times (s) — reasonable estimates if no mat file exists
if use_hover
    FALLBACK = [12; 17; 30; 35];   % [t1, t2, t3, t4] end times
else
    FALLBACK = [12; 28; 33];       % [t1, t3, t4] end times (no hover)
end

% --- No path supplied ---
if nargin < 1 || isempty(ref_traj_path)
    warning('load_phase_times: no path given — using fallback times');
    t_phase_nominal = FALLBACK;
    return;
end

% --- File missing ---
if ~isfile(ref_traj_path)
    warning('load_phase_times: file not found: %s\n  Using fallback.', ref_traj_path);
    t_phase_nominal = FALLBACK;
    return;
end

% --- Load ---
try
    data = load(ref_traj_path);
catch ME
    warning('load_phase_times: load failed: %s', ME.message);
    t_phase_nominal = FALLBACK;
    return;
end

% --- Extract from phase_idx (preferred) ---
if isfield(data, 'phase_idx') && isfield(data, 't_ref')
    t_ref     = data.t_ref(:);
    phase_idx = data.phase_idx(:);

    % Active phases depend on hover flag
    if use_hover
        active_phases = [1 2 3 4];
    else
        active_phases = [1 3 4];
    end

    t_ends = zeros(numel(active_phases), 1);
    for k = 1:numel(active_phases)
        ph  = active_phases(k);
        idx = find(phase_idx == ph, 1, 'last');
        if ~isempty(idx)
            t_ends(k) = t_ref(idx);
        else
            warning('load_phase_times: phase %d not found — using fallback', ph);
            t_ends(k) = FALLBACK(k);
        end
    end
    t_phase_nominal = t_ends;

% --- Fallback: detect from altitude profile ---
elseif isfield(data, 't_ref') && isfield(data, 'x_ref')
    warning('load_phase_times: phase_idx missing — detecting from altitude');
    t_ref = data.t_ref(:);
    Z_ref = data.x_ref(:, 3);

    % t1_end: Z first reaches 295 m
    idx = find(Z_ref >= 295, 1, 'first');
    t1_end = t_ref(max(1, idx));

    if use_hover
        t2_end = t1_end + 5.0;
        idx_post = find(t_ref > t2_end, 1, 'first');
    else
        idx_post = find(t_ref > t1_end, 1, 'first');
    end

    % t3_end: Z drops below 8 m after hover/ascent
    idx34 = find(Z_ref(idx_post:end) < 8, 1, 'first');
    t3_end = t_ref(idx_post + idx34 - 1);
    t4_end = t3_end + 5.0;

    if use_hover
        t_phase_nominal = [t1_end; t2_end; t3_end; t4_end];
    else
        t_phase_nominal = [t1_end; t3_end; t4_end];
    end

else
    warning('load_phase_times: unrecognised mat file — using fallback');
    t_phase_nominal = FALLBACK;
end

fprintf('Phase nominal end times: %s\n', mat2str(t_phase_nominal', 4));
end