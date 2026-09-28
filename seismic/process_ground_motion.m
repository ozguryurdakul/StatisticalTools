function R = process_ground_motion(time, acc_in, opt)
%PROCESS_GROUND_MOTION  Baseline-correct, filter and integrate an accelerogram.
%
%   R = process_ground_motion(time, acc, opt)
%
%   Processing chain (every step acts on ACCELERATION only, so velocity and
%   displacement are always the exact integrals of the acceleration that is
%   reported -- a -> v -> d stay mutually consistent):
%
%     1. unit conversion to cm/s^2
%     2. mean removal (pre-event mean if opt.pre_event_s > 0, else whole record)
%     3. least-squares polynomial baseline on acceleration (order 0..3,
%        same idea as SeismoSignal's "baseline correction")
%     4. cosine taper at both ends
%     5. zero padding (Boore 2005 / Boore & Bommer 2005 length)
%     6. Butterworth filter in second-order sections (numerically stable),
%        acausal (forward-backward) or causal
%     7. integration with initial conditions v(0)=0, d(0)=0 at the start of
%        the PADDED record; the pads are removed only after integration
%     8. optional displacement-polynomial correction (applied back to acc)
%     9. PGA / PGV / PGD, Fourier amplitude spectrum, drift diagnostics
%
%   Output units: acceleration cm/s^2, velocity cm/s, displacement cm.
%
%   opt fields (all optional; defaults in brackets):
%     units          'g' | 'cm/s2' | 'm/s2' | 'mm/s2'              ['cm/s2']
%     pre_event_s    length of pre-event noise used for the mean [s]  [0]
%     baseline_order polynomial order fitted to acc, -1 = none        [1]
%     taper_s        cosine-taper length at each end [s]              [1]
%     filter_type    'highpass' | 'bandpass' | 'none'                 ['bandpass']
%     fc_low         high-pass corner [Hz]                            [0.1]
%     fc_high        low-pass corner  [Hz] (bandpass only)            [25]
%     order          Butterworth order PER PASS                       [4]
%     causal         false = zero-phase forward-backward              [false]
%     pad            true = Boore zero pads, false = none             [true]
%     keep_pads      true = return the padded series too              [false]
%     integration    'linear' (exact for piecewise-linear acc) | 'trapz' ['linear']
%     disp_poly_order  0 = off; else fit t^2..t^N to displacement and
%                      remove its 2nd derivative from acc (PEER-style)  [0]

if nargin < 3, opt = struct(); end
opt = set_defaults(opt);

time = time(:); acc = acc_in(:);
n    = numel(acc);
if numel(time) ~= n, error('time and acceleration must have equal length'); end

% --- sampling interval: take it from the file, and make sure it is uniform
dts = diff(time);
dt  = median(dts);
if max(abs(dts - dt)) > 1e-3*dt
    error('Time vector is not uniformly sampled (dt varies by %.3g s). Resample first.', ...
          max(abs(dts - dt)));
end
Fs = 1/dt;

% --- 1. units -> cm/s^2
switch lower(opt.units)
    case 'g',      k = 980.665;
    case 'cm/s2',  k = 1;
    case 'm/s2',   k = 100;
    case 'mm/s2',  k = 0.1;
    otherwise, error('Unknown units "%s"', opt.units);
end
acc = acc * k;
acc_raw = acc;

% --- 2. mean removal
npre = round(opt.pre_event_s / dt);
if npre >= 10
    acc = acc - mean(acc(1:min(npre, n)));
else
    acc = acc - mean(acc);
end

% --- 3. polynomial baseline on acceleration (least squares)
%     time is centred and scaled so that polyfit stays well conditioned
if opt.baseline_order >= 0
    tt = (0:n-1)' * dt;
    x  = (tt - mean(tt)) / std(tt);
    p  = polyfit(x, acc, opt.baseline_order);
    acc = acc - polyval(p, x);
end

% --- 4. cosine taper
m = min(round(opt.taper_s / dt), floor(n/2));
if m > 1
    w = 0.5 * (1 - cos(pi * (0:m-1)' / m));
    acc(1:m)         = acc(1:m)         .* w;
    acc(end-m+1:end) = acc(end-m+1:end) .* flipud(w);
end

% --- 5. zero pads (Boore 2005): total pad = 1.5*order/fc for an acausal
%     filter, split evenly between the two ends.
use_filter = ~strcmpi(opt.filter_type, 'none');
npad = 0;
if use_filter && opt.pad
    npad = ceil(0.5 * 1.5 * opt.order / opt.fc_low / dt);
end
acc_p = [zeros(npad,1); acc; zeros(npad,1)];

% --- 6. Butterworth filter in SOS form
%     [b,a] form loses accuracy for low Wn and is UNSTABLE for order >= 6
%     at fc ~ 0.1 Hz / Fs = 200 Hz, so zpk -> sos is used instead.
if use_filter
    nyq = Fs/2;
    switch lower(opt.filter_type)
        case 'highpass'
            [z, pz, g] = butter(opt.order, opt.fc_low/nyq, 'high');
        case 'bandpass'
            if opt.fc_high >= nyq
                error('fc_high (%.3g Hz) must be below Nyquist (%.3g Hz)', opt.fc_high, nyq);
            end
            [z, pz, g] = butter(opt.order, [opt.fc_low opt.fc_high]/nyq);
        otherwise
            error('Unknown filter_type "%s"', opt.filter_type);
    end
    sos = zp2sos(z, pz, g);
    acc_p = sosfilt(sos, acc_p);
    if ~opt.causal
        % forward-backward pass: zero phase, magnitude |H|^2
        % (effective order = 2*opt.order, -6 dB at the corner)
        acc_p = flipud(sosfilt(sos, flipud(acc_p)));
    end
end

% --- 7. integrate the whole padded record, starting from rest
[vel_p, dis_p] = integrate_acc(acc_p, dt, opt.integration);

% --- 8. optional displacement-polynomial correction (kept consistent by
%     subtracting the 2nd derivative of the fitted polynomial from acc)
if opt.disp_poly_order >= 2
    np_ = numel(acc_p);
    tp  = (0:np_-1)' * dt;
    T   = tp(end);
    s   = tp / T;                         % scaled time 0..1
    pw  = 2:opt.disp_poly_order;
    A   = s .^ pw;                        % no constant / linear term
    c   = A \ dis_p;
    % d2/dt2 of c_k (t/T)^k = c_k k (k-1) t^(k-2) / T^k
    a_corr = (s .^ (pw-2)) * (c(:) .* (pw(:).*(pw(:)-1))) / T^2;
    acc_p  = acc_p - a_corr;
    [vel_p, dis_p] = integrate_acc(acc_p, dt, opt.integration);
end

% --- remove pads (after integration!)
idx = (npad+1):(npad+n);
R.time = time;
R.dt   = dt;
R.acc_raw = acc_raw;
R.acc  = acc_p(idx);
R.vel  = vel_p(idx);
R.dis  = dis_p(idx);
R.npad = npad;
if opt.keep_pads
    R.time_padded = time(1) + ((0:numel(acc_p)-1)' - npad) * dt;
    R.acc_padded  = acc_p;
    R.vel_padded  = vel_p;
    R.dis_padded  = dis_p;
end

% --- 9. peak values
[R.PGA, i] = max(abs(R.acc)); R.t_PGA = time(i);
[R.PGV, i] = max(abs(R.vel)); R.t_PGV = time(i);
[R.PGD, i] = max(abs(R.dis)); R.t_PGD = time(i);

% Fourier amplitude spectrum, |FFT|*dt  [cm/s], one-sided, odd/even-safe
L = numel(R.acc);
Y = fft(R.acc);
nh = floor(L/2) + 1;
R.freq = (0:nh-1)' * Fs / L;
R.FAS  = abs(Y(1:nh)) * dt;

% drift diagnostics: after a good correction the tail velocity ~ 0 and the
% tail displacement is flat (no ramp). Uses the last 10 % of the record.
nt = max(10, round(0.1*n));
R.diag.v_end         = R.vel(end);
R.diag.v_tail_mean   = mean(R.vel(end-nt+1:end));
pd = polyfit((0:nt-1)'*dt, R.dis(end-nt+1:end), 1);
R.diag.d_tail_slope  = pd(1);            % cm/s, should be ~0
R.diag.d_end         = R.dis(end);
R.diag.v_start       = R.vel(1);         % non-zero only because of acausal pads
R.diag.d_start       = R.dis(1);
R.opt = opt;
end

% -------------------------------------------------------------------------
function [v, d] = integrate_acc(a, dt, method)
a = a(:);
v = [0; cumsum((a(1:end-1) + a(2:end)) * dt/2)];
switch lower(method)
    case 'linear'
        % exact double integral of a piecewise-linear acceleration
        dd = v(1:end-1)*dt + (2*a(1:end-1) + a(2:end)) * dt^2/6;
        d  = [0; cumsum(dd)];
    case 'trapz'
        d  = [0; cumsum((v(1:end-1) + v(2:end)) * dt/2)];
    otherwise
        error('Unknown integration "%s"', method);
end
end

function opt = set_defaults(opt)
def = struct('units','cm/s2', 'pre_event_s',0, 'baseline_order',1, ...
             'taper_s',1, 'filter_type','bandpass', 'fc_low',0.1, ...
             'fc_high',25, 'order',4, 'causal',false, 'pad',true, ...
             'keep_pads',false, 'integration','linear', 'disp_poly_order',0);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opt, f{i}), opt.(f{i}) = def.(f{i}); end
end
end
