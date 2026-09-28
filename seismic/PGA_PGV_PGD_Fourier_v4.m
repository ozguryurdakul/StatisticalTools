clear
close all
clc
%% PGA / PGV / PGD and Fourier spectrum -- v4 (audited, see AUDIT_PGV_PGD.md)
%
% Input: 2-column text file [time(s)  acceleration].
% All processing is in process_ground_motion.m. Output units are
% cm/s^2, cm/s and cm (NOT mm).

%% ---- Settings ------------------------------------------------------------
file = 'kocaeli_ew.txt';

opt.units          = 'cm/s2';   % units of the acceleration column: 'g','cm/s2','m/s2'
opt.pre_event_s    = 0;         % >0: remove mean of first N s instead of whole-record mean
opt.baseline_order = 1;         % polynomial fitted to acc: 0 const, 1 linear, 2, 3; -1 none
opt.taper_s        = 1;         % cosine taper at each end [s]
opt.filter_type    = 'bandpass';% 'highpass' | 'bandpass' | 'none'
opt.fc_low         = 0.1;       % [Hz] PGD is very sensitive to this - see audit
opt.fc_high        = 25;        % [Hz]
opt.order          = 4;         % Butterworth order per pass
opt.causal         = false;     % false = zero-phase (acausal) with Boore zero pads
opt.integration    = 'linear';  % 'linear' (exact for piecewise-linear acc) or 'trapz'
opt.disp_poly_order = 0;        % e.g. 6 to remove residual displacement drift; 0 = off

%% ---- Read ----------------------------------------------------------------
data = importdata(file);
if isstruct(data), data = data.data; end
time = data(:,1);
acc  = data(:,2);

%% ---- Uncorrected reference (to see what the corrections do) ---------------
ref  = opt; ref.baseline_order = -1; ref.filter_type = 'none'; ref.taper_s = 0;
R0   = process_ground_motion(time, acc, ref);    % mean removal only
R    = process_ground_motion(time, acc, opt);

fprintf('dt = %.4f s  (Fs = %.1f Hz), N = %d, zero pads = %d samples per end\n', ...
        R.dt, 1/R.dt, numel(R.time), R.npad);
fprintf('PGA = %9.3f cm/s^2  (%.4f g)  at t = %.3f s\n', R.PGA, R.PGA/980.665, R.t_PGA);
fprintf('PGV = %9.3f cm/s             at t = %.3f s\n', R.PGV, R.t_PGV);
fprintf('PGD = %9.3f cm               at t = %.3f s\n', R.PGD, R.t_PGD);
fprintf('Drift checks: v_end = %.3f cm/s, tail displacement slope = %.4f cm/s\n', ...
        R.diag.v_end, R.diag.d_tail_slope);

%% ---- PGD sensitivity to the high-pass corner (always report this) ---------
fcs = [0.03 0.05 0.07 0.1 0.15 0.2];
fprintf('\n fc_low [Hz]   PGV [cm/s]   PGD [cm]\n');
for f = fcs
    o = opt; o.fc_low = f;
    Rs = process_ground_motion(time, acc, o);
    fprintf('   %5.3f      %8.2f    %8.2f\n', f, Rs.PGV, Rs.PGD);
end

%% ---- Plots ----------------------------------------------------------------
figure(1)
subplot(3,1,1)
plot(R.time, R.acc, 'k'); grid on
ylabel('Acc. [cm/s^2]'); title(sprintf('Corrected, PGA = %.2f cm/s^2', R.PGA))
subplot(3,1,2)
plot(R0.time, R0.vel, 'color', [0.7 0.7 0.7]); hold on
plot(R.time, R.vel, 'k'); grid on
ylabel('Vel. [cm/s]'); legend('mean removal only', 'corrected', 'Location', 'best')
title(sprintf('PGV = %.2f cm/s', R.PGV))
subplot(3,1,3)
plot(R0.time, R0.dis, 'color', [0.7 0.7 0.7]); hold on
plot(R.time, R.dis, 'k'); grid on
xlabel('Time [s]'); ylabel('Disp. [cm]')
title(sprintf('PGD = %.2f cm', R.PGD))

figure(2)
loglog(R0.freq(2:end), R0.FAS(2:end), 'color', [0.7 0.7 0.7]); hold on
loglog(R.freq(2:end), R.FAS(2:end), 'k'); grid on
xlabel('Frequency [Hz]'); ylabel('Fourier amplitude [cm/s]')
legend('mean removal only', 'corrected', 'Location', 'best')
title('Fourier amplitude spectrum of acceleration')
xlim([0.01 1/(2*R.dt)])

%% ---- Export (with header) --------------------------------------------------
fid = fopen('results_corrected.txt', 'w');
fprintf(fid, '%% PGA=%.6g cm/s2  PGV=%.6g cm/s  PGD=%.6g cm\n', R.PGA, R.PGV, R.PGD);
fprintf(fid, '%% time[s]\tacc[cm/s2]\tvel[cm/s]\tdisp[cm]\n');
fprintf(fid, '%.5f\t%.6e\t%.6e\t%.6e\n', [R.time R.acc R.vel R.dis]');
fclose(fid);
