clear all
close all
clc

%% Input data (using importdata)

data=importdata('kocaeli_ew.txt');
T= 0.005            % Sampling period 
Fs = 1/T;                       % Sampling frequency [Hz]

time=data(:,1);
acceleration=data(:,2); %units in  cm2/s

dt = mean(diff(time));

%% Data preprocessing
% Remove any linear trend from acceleration
acc_detrended = detrend(acceleration);

% Optional: High-pass filter to remove low-frequency noise/drift
fc = 0.1; % Cutoff frequency [Hz] - adjust based on your signal
if fc < Fs/2 && fc > 0
    [b,a] = butter(4, fc/(Fs/2), 'high');
    acc_filtered = filtfilt(b, a, acc_detrended);
else
    acc_filtered = acc_detrended;
end

%% --- Method 1: Simple cumtrapz double integration ---
velocity1 = cumtrapz(time, acceleration);
displacement1 = cumtrapz(time, velocity1);

%% --- Method 2: Double integration with zero-mean acceleration ---
acc_zero_mean = acceleration - mean(acceleration);
velocity2 = cumtrapz(time, acc_zero_mean);
displacement2 = cumtrapz(time, velocity2);

%% --- Method 3: Improved method with detrended acceleration ---
velocity3 = cumtrapz(time, acc_detrended);
displacement3 = cumtrapz(time, velocity3);

%% --- Method 4: Best method - filtered + detrended + baseline correction ---
velocity4 = cumtrapz(time, acc_filtered);
velocity4 = detrend(velocity4); % Remove linear drift from velocity
displacement4 = cumtrapz(time, velocity4);
displacement4 = detrend(displacement4); % Remove linear drift from displacement

%% Convert for plotting
% Method 1
velocity1_cm = velocity1; % cm/s
displacement1_mm = displacement1 ; % mm

% Method 2
velocity2_cm = velocity2;
displacement2_mm = displacement2 ;

% Method 3
velocity3_cm = velocity3;
displacement3_mm = displacement3 ;

% Method 4 (improved)
velocity4_cm = velocity4;
displacement4_mm = displacement4 ;


%% --- Plots ---
% Acceleration comparison
figure(1)
subplot(2,1,1)
plot(time, acceleration, 'b', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
title('Original Acceleration vs Time');
grid on

subplot(2,1,2)
plot(time, acc_filtered, 'r', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
title('Filtered & Detrended Acceleration vs Time');
grid on

% Velocity comparison
figure(2)
plot(time, velocity1_cm, 'b', ...
     time, velocity2_cm, 'g-.', ...
     time, velocity3_cm, 'r--', ...
     time, velocity4_cm, 'k', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Velocity [cm/s]');
legend('Method 1: Simple', 'Method 2: Zero-mean', ...
       'Method 3: Detrended', 'Method 4: Improved', 'Location', 'best');
title('Velocity vs Time - Method Comparison');
grid on

% Displacement comparison
figure(3)
plot(time, displacement1_mm, 'b', ...
     time, displacement2_mm, 'g-.', ...
     time, displacement3_mm, 'r--', ...
     time, displacement4_mm, 'k', 'LineWidth', 1.2);
xlabel('Time [s]');
ylabel('Displacement [mm]');
legend('Method 1: Simple', 'Method 2: Zero-mean', ...
       'Method 3: Detrended', 'Method 4: Improved', 'Location', 'best');
title('Displacement vs Time - Method Comparison');
grid on

%% --- Fourier transformation of acceleration ---
figure(4)
L = length(acceleration);
Y_orig = fft(acceleration);
Y_filt = fft(acc_filtered);

P2_orig = abs(Y_orig/L);
P1_orig = P2_orig(1:L/2+1);
P1_orig(2:end-1) = 2*P1_orig(2:end-1);

P2_filt = abs(Y_filt/L);
P1_filt = P2_filt(1:L/2+1);
P1_filt(2:end-1) = 2*P1_filt(2:end-1);

f = Fs*(0:(L/2))/L;

subplot(2,1,1)
plot(f, P1_orig, 'b', 'LineWidth', 1.2);
xlabel('Frequency [Hz]');
ylabel('|Amplitude (m/s²)|');
title('Original Acceleration - Frequency Spectrum');
grid on
xlim([0 min(50, Fs/2)]); % Limit to 50 Hz or Nyquist frequency

subplot(2,1,2)
plot(f, P1_filt, 'r', 'LineWidth', 1.2);
xlabel('Frequency [Hz]');
ylabel('|Amplitude (m/s²)|');
title('Filtered Acceleration - Frequency Spectrum');
grid on
xlim([0 min(50, Fs/2)]); % Limit to 50 Hz or Nyquist frequency

%% --- Statistical analysis ---
figure(5)
subplot(2,2,1)
histogram(acceleration, 30, 'Normalization', 'probability');
title('Original Acceleration Distribution');
xlabel('Acceleration [m/s²]');
ylabel('Probability');
grid on

subplot(2,2,2)
histogram(velocity4_cm, 30, 'Normalization', 'probability');
title('Velocity Distribution (Method 4)');
xlabel('Velocity [cm/s]');
ylabel('Probability');
grid on

subplot(2,2,3)
histogram(displacement4_mm, 30, 'Normalization', 'probability');
title('Displacement Distribution (Method 4)');
xlabel('Displacement [mm]');
ylabel('Probability');
grid on

subplot(2,2,4)
% Time series of all three signals
yyaxis left
plot(time, acceleration, 'b', 'LineWidth', 0.8);
ylabel('Acceleration [m/s²]');
yyaxis right
plot(time, displacement4_mm, 'r', 'LineWidth', 1);
ylabel('Displacement [mm]');
xlabel('Time [s]');
title('Acceleration vs Displacement');
grid on

%% --- Export results ---
% Export all methods for comparison
results_all = [time, displacement1_mm, displacement2_mm, displacement3_mm, displacement4_mm];
results_header = 'Time[s], Disp1[mm], Disp2[mm], Disp3[mm], Disp4[mm]';

% Save comprehensive results
save('results_all_methods.txt', 'results_all', '-ascii');

% Save best method only
results_best = [time, velocity4_cm, displacement4_mm];
save('results_best_method.txt', 'results_best', '-ascii');

% Analysis complete! Method 4 is recommended for most accurate results.