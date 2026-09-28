function test_process_ground_motion()
%TEST_PROCESS_GROUND_MOTION  Checks against records with an analytic answer.
%   Run:  test_process_ground_motion   (Octave: pkg load signal first)

dt = 0.005; t = (0:dt:60-dt)';
[a, v, d] = pulse(t, 20, 3, 0.3, 30);          % D = 30 cm, 0.3 Hz pulse
nfail = 0;

% 1. Integration alone: clean record, no processing -> must be exact
o = struct('baseline_order', -1, 'filter_type', 'none', 'taper_s', 0);
R = process_ground_motion(t, a - 0, o);
% process removes the mean; the analytic a has ~zero mean already
nfail = nfail + check('integration PGV', R.PGV, max(abs(v)), 1e-3);
nfail = nfail + check('integration PGD', R.PGD, max(abs(d)), 1e-3);

% 2. Units: same record given in g must give identical cm results
Rg = process_ground_motion(t, a/980.665, setfield(o, 'units', 'g'));
nfail = nfail + check('units g -> cm PGD', Rg.PGD, R.PGD, 1e-9);

% 3. Baseline errors (offset + step + noise), default processing
randn('state', 0); %#ok<RAND>
ac = a + 0.3 + 0.8*(t > 18) + 0.3*randn(size(t));
Rc = process_ground_motion(t, ac, struct());
nfail = nfail + check('corrected PGV', Rc.PGV, max(abs(v)), 0.02);
nfail = nfail + check('corrected PGD', Rc.PGD, max(abs(d)), 0.03);

% 4. Internal consistency: d is the integral of v, v the integral of a
dv = diff(Rc.dis)/dt - (Rc.vel(1:end-1) + Rc.vel(2:end))/2;
nfail = nfail + check('d/v consistency', max(abs(dv)), 0, 0, 1e-3);

% 5. High-order filter (would blow up with [b,a] form)
Rh = process_ground_motion(t, ac, struct('order', 8, 'filter_type', 'highpass'));
nfail = nfail + check('order-8 highpass PGD', Rh.PGD, max(abs(d)), 0.05);

% 6. Odd-length record (v3 FFT indexing crashed here)
Ro = process_ground_motion(t(1:end-1), ac(1:end-1), struct());
nfail = nfail + check('odd-length FAS size', numel(Ro.FAS), floor(numel(t(1:end-1))/2)+1, 0);

if nfail == 0, fprintf('ALL TESTS PASSED\n'); else, error('%d test(s) failed', nfail); end
end

function f = check(name, got, want, rtol, atol)
if nargin < 5, atol = 0; end
ok = abs(got - want) <= rtol*abs(want) + atol;
fprintf('%-24s got %12.5g  expected %12.5g  %s\n', name, got, want, ifelse(ok));
f = ~ok;
end

function s = ifelse(ok)
if ok, s = 'ok'; else, s = 'FAIL'; end
end

function [a, v, d] = pulse(t, t0, s, f, D)
g = exp(-((t-t0)/s).^2); w = 2*pi*f;
gp  = -2*(t-t0)/s^2 .* g;
gpp = (4*(t-t0).^2/s^4 - 2/s^2) .* g;
d = D * g .* sin(w*t);
v = D * (gp .* sin(w*t) + g*w .* cos(w*t));
a = D * (gpp .* sin(w*t) + 2*gp*w .* cos(w*t) - g*w^2 .* sin(w*t));
end
