function out = run_nominal(stop_time)
% RUN_NOMINAL  Single entry point for the closed-loop sim.
% L. White, University of Sheffield May 2026
%   out = run_nominal()      runs for the full reference duration
%   out = run_nominal(10)    runs for 10 seconds
if nargin < 1, stop_time = []; end


% 1. Initialise workspace
params    = plant_params();
assert(isfield(params,'reference') && ~isempty(params.reference), ...
    'plant_params did not load reference_trajectory.mat. Fix the path inside plant_params.m.');


reference = params.reference;
assignin('base','params',    params);
assignin('base','reference', reference);

% 2. Reset every persistent variable in every wrapper
clear ekf_wrapper scheduler_wrapper guidance_wrapper attitude_wrapper actuator_wrapper sensors_wrapper
clear plant_dynamics imu_model gps_model baro_model wind_model

% 3. Pick stop time
if isempty(stop_time)
    stop_time = reference.t(end) + 1;     % full mission + 1 s margin
end
fprintf('Running vtvl_sim for %.2f s ...\n', stop_time);

% 4. Run the sim, capture full error info
try
    out = sim('vtvl_sim', ...
        'StopTime',       num2str(stop_time), ...
        'SimulationMode', 'normal', ...
        'ReturnWorkspaceOutputs','on');
    fprintf('Sim completed: reached t=%.3f s (requested %.2f s)\n', ...
        out.tout(end), stop_time);

% 5. Pull all logged signals 
flds = fieldnames(out);
fprintf('\nFields in out:  '); fprintf('%s  ', flds{:}); fprintf('\n');

xt = squeeze(out.x_true_log); if size(xt,1)<size(xt,2), xt = xt'; end
xh = squeeze(out.x_hat_log);  if size(xh,1)<size(xh,2), xh = xh'; end
Tc = squeeze(out.T_cmd_log);  Tc = Tc(:);
al = squeeze(out.alpha_log);  al = al(:);
ph = squeeze(out.phase_log);  ph = ph(:);

% Matched time vectors (each block has its own sample rate)
N_xt = size(xt,1); t_xt = linspace(0, out.tout(end), N_xt)';
N_xh = size(xh,1); t_xh = linspace(0, out.tout(end), N_xh)';
N_Tc = numel(Tc);  t_Tc = linspace(0, out.tout(end), N_Tc )';
N_al = numel(al);  t_al = linspace(0, out.tout(end), N_al )';
N_ph = numel(ph);  t_ph = linspace(0, out.tout(end), N_ph )';

fprintf('Sample counts: x_true=%d  x_hat=%d  T_cmd=%d  alpha=%d  phase=%d\n',...
        N_xt, N_xh, N_Tc, N_al, N_ph);

% Convert true and MEKF quaternions to Euler (deg)
quat2eul_deg = @(Q) [ ...
    asind(max(-1,min(1, 2*(Q(:,1).*Q(:,3) - Q(:,4).*Q(:,2))))), ...     % theta
    atan2d(2*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)), 1-2*(Q(:,2).^2+Q(:,3).^2))]; % phi
eul_true = quat2eul_deg(xt(:,8:11));    % cols: theta, phi (deg)
eul_hat  = quat2eul_deg(xh(:,7:10));

% Compare MEKF attitude to truth
N_min = min(size(xt,1), size(xh,1));
fprintf('Max |theta_true - theta_hat|: %.3f deg\n', ...
    max(abs(eul_true(1:N_min,1) - eul_hat(1:N_min,1))));
fprintf('Max |phi_true   - phi_hat|:   %.3f deg\n', ...
    max(abs(eul_true(1:N_min,2) - eul_hat(1:N_min,2))));

% Reference body Euler for comparison
eul_ref = quat2eul_deg(reference.x(:,8:11));

% 6. Summary numbers
fprintf('\n--- Summary ---\n');
fprintf('Final pos truth : [%+8.2f %+8.2f %+8.2f] m  (ref pad %+d %+d 0)\n', ...
        xt(end,1:3), params.r_pad(1), params.r_pad(2));
fprintf('Max altitude    : %.2f m  (ref apex 300)\n', max(xt(:,3)));
fprintf('Final mass      : %.3f kg  (m0=%.2f, ms=%.2f)\n', xt(end,7), params.m0, params.ms);
fprintf('Time mass hits ms: ');
i_dry = find(xt(:,7) <= params.ms*1.001, 1, 'first');
if isempty(i_dry), fprintf('never\n');
else, fprintf('%.2f s\n', t_xt(i_dry)); end
fprintf('Max |theta_true|: %.2f deg,  max |phi_true|: %.2f deg\n', ...
        max(abs(eul_true(:,1))), max(abs(eul_true(:,2))));
fprintf('Alpha sat (>=99%% limit): %.1f%% of samples\n', ...
        100*mean(abs(al) >= params.alpha_max*0.99));

% 7. Plots

% Note: x_hat layout = [r(3); v(3); q(4); b_g(3); b_a(3); m(1)]
idx_0   = 1;
idx_500 = find(t_xt >= 0.5, 1, 'first');
fprintf('---- EKF snapshot ----\n');
fprintf('t=0:    q_hat=[% .4f % .4f % .4f % .4f]  b_g=[% .4f % .4f % .4f]  b_a=[% .4f % .4f % .4f]\n', ...
    xh(idx_0,   7:10), xh(idx_0,   11:13), xh(idx_0,   14:16));
fprintf('t=0.5:  q_hat=[% .4f % .4f % .4f % .4f]  b_g=[% .4f % .4f % .4f]  b_a=[% .4f % .4f % .4f]\n', ...
    xh(idx_500, 7:10), xh(idx_500, 11:13), xh(idx_500, 14:16));
fprintf('        q_true=[% .4f % .4f % .4f % .4f] (for comparison)\n', xt(idx_500, 8:11));

% Proper unwrapped diff in deg
d_theta = abs(rad2deg(wrapToPi(deg2rad(eul_true(:,1) - eul_hat(:,1)))));
d_phi   = abs(rad2deg(wrapToPi(deg2rad(eul_true(:,2) - eul_hat(:,2)))));
fprintf('Unwrapped max |theta_true - theta_hat|: %.3f deg\n', max(d_theta));
fprintf('Unwrapped max |phi_true   - phi_hat|:   %.3f deg\n', max(d_phi));

figure('Name','Trajectory','NumberTitle','off','Position',[100 100 1000 700]);
subplot(2,2,1); plot3(xt(:,1),xt(:,2),xt(:,3),'b'); hold on; grid on;
    plot3(reference.x(:,1),reference.x(:,2),reference.x(:,3),'r--');
    plot3(params.r_pad(1),params.r_pad(2),0,'kx','MarkerSize',12,'LineWidth',2);
    xlabel('X');ylabel('Y');zlabel('Z'); title('3-D path'); axis equal;
    legend('Sim','Ref','Pad');
subplot(2,2,2); plot(t_xt,xt(:,3),'b', reference.t,reference.x(:,3),'r--');
    grid on; ylabel('Z (m)'); xlabel('t'); title('Altitude');
subplot(2,2,3); plot(t_xt,xt(:,1),'b', reference.t,reference.x(:,1),'r--');
    grid on; ylabel('X (m)'); xlabel('t'); title('X position');
subplot(2,2,4); plot(t_xt,xt(:,2),'b', reference.t,reference.x(:,2),'r--');
    grid on; ylabel('Y (m)'); xlabel('t'); title('Y position');

figure('Name','Body attitude','NumberTitle','off','Position',[100 100 1000 500]);
subplot(1,2,1); plot(t_xt,eul_true(:,1),'b', reference.t,eul_ref(:,1),'r--');
    grid on; ylabel('\theta (deg)'); xlabel('t'); title('Pitch');
    legend('Sim truth','Ref','Location','best');
    yline(30,'k:'); yline(-30,'k:');
subplot(1,2,2); plot(t_xt,eul_true(:,2),'b', reference.t,eul_ref(:,2),'r--');
    grid on; ylabel('\phi (deg)'); xlabel('t'); title('Roll');
    yline(30,'k:'); yline(-30,'k:');

figure('Name','Actuators and mass','NumberTitle','off','Position',[100 100 1000 700]);
subplot(2,2,1); plot(t_Tc, Tc,'b', reference.t,reference.u(:,1),'r--');
    grid on; ylabel('T cmd (N)'); xlabel('t'); title('Thrust command');
    yline(params.Tmin,'k:'); yline(params.Tmax,'k:');
subplot(2,2,2); plot(t_al, al*180/pi,'b', reference.t,reference.u(:,2),'r--');
    grid on; ylabel('\alpha (deg)'); xlabel('t'); title('Gimbal pitch (cmd)');
    yline(params.alpha_max*180/pi,'k:'); yline(-params.alpha_max*180/pi,'k:');
subplot(2,2,3); plot(t_xt, xt(:,7),'b', reference.t,reference.x(:,7),'r--');
    grid on; ylabel('mass (kg)'); xlabel('t'); title('Vehicle mass');
    yline(params.ms,'k:','dry');
subplot(2,2,4); plot(t_ph, ph,'b', reference.t, double(reference.phase),'r--');
    grid on; ylabel('phase'); xlabel('t'); title('Phase'); ylim([0 5]);
catch ME
    fprintf('\nCouldn''t auto-extract x_true_log: %s\n', ME.message);
    fprintf('Tell me the actual signal names in your model.\n');
end
end

