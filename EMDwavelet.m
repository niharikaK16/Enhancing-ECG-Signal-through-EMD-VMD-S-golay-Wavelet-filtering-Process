% Adding Gaussian noise to ECG signal
function noisy_ecg = add_gaussian_noise(Original_signal, snr_db)
signal_power = rms(Original_signal)^2;
noise_power = signal_power / (10^(snr_db/10));
noise = sqrt(noise_power) * randn(size(Original_signal));
noisy_ecg = Original_signal + noise;
end
 
% Load ECG signal data
load('/MATLAB Drive/100m.mat');
 
val = (val - 0) / 200; 
 
sig = val(1, 1:3600); 
 
Fs = 360; 
 
t = (0:length(sig)-1) / Fs; 
 
snr_db = 10; 
noisy_sig = add_gaussian_noise(sig, snr_db);
 
% Plot original and noisy ECG signals
figure;
subplot(2,1,1);
plot(t, sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Original ECG Signal');
grid on;
 
subplot(2,1,2);
plot(t, noisy_sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Noisy ECG Signal');
grid on;


% Empirical Mode Decomposition (EMD) to decompose the signal into IMFs
[IMFs, residue] = emd(noisy_sig);
 
% Denoise each IMF using Wavelet thresholding
numIMFs = size(IMFs, 2);
denoisedIMFs = zeros(size(IMFs));
 
for k = 1:numIMFs
   [c, l] = wavedec(IMFs(:, k), 5, 'db4'); 
   sigma = median(abs(c)) / 0.6745; 
   thr = sigma * sqrt(2 * log(length(c))); 
   c = wthresh(c, 's', thr); 
   denoisedIMFs(:, k) = waverec(c, l, 'db4');
end
 
% Reconstruct the denoised signal from the denoised IMFs
denoised_sig = sum(denoisedIMFs, 2) + residue;
 
% Plot denoised ECG signal
figure;
plot(t, denoised_sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Denoised ECG Signal (EMD + Wavelet)');
grid on;